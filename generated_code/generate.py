#!/usr/bin/env python3
##
# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2019, SeisSol Group
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
#

import sys
import argparse
import importlib.util
import inspect

from yateto import useArchitectureIdentifiedBy, Generator
from yateto import gemm_configuration
from yateto.gemm_configuration import GeneratorCollection, LIBXSMM, PSpaMM, GemmForge
from yateto.gemm_configuration import GeneratorCollection, MKL, BLIS, OpenBLAS

import DynamicRupture
import Plasticity
import SurfaceDisplacement
import Point
import memlayout

cmdLineParser = argparse.ArgumentParser()
cmdLineParser.add_argument('--equations')
cmdLineParser.add_argument('--matricesDir')
cmdLineParser.add_argument('--outputDir')
cmdLineParser.add_argument('--host_arch')
cmdLineParser.add_argument('--compute_arch', default=None)
cmdLineParser.add_argument('--order', type=int)
cmdLineParser.add_argument('--numberOfMechanisms', type=int)
cmdLineParser.add_argument('--memLayout')
cmdLineParser.add_argument('--multipleSimulations', type=int)
cmdLineParser.add_argument('--dynamicRuptureMethod')
cmdLineParser.add_argument('--PlasticityMethod')
cmdLineParser.add_argument('--gemm_tools')
cmdLineArgs = cmdLineParser.parse_args()

if cmdLineArgs.memLayout == 'auto':
  # TODO(Lukas) Don't hardcode this
  env = {
    'equations': cmdLineArgs.equations,
    'order': cmdLineArgs.order,
    'arch': cmdLineArgs.arch,
    'multipleSimulations': cmdLineArgs.multipleSimulations
  }
  mem_layout = memlayout.guessMemoryLayout(env)
else:
  mem_layout = cmdLineArgs.memLayout
  

arch = useArchitectureIdentifiedBy(cmdLineArgs.host_arch, cmdLineArgs.compute_arch)

equationsSpec = importlib.util.find_spec(cmdLineArgs.equations)
try:
  equations = equationsSpec.loader.load_module()
except:
  raise RuntimeError('Could not find kernels for ' + cmdLineArgs.equations)

gpu_platforms = ['nvidia']
platforms = ['gpu', 'cpu'] if cmdLineArgs.compute_arch[1:] in gpu_platforms else ['cpu']

adgArgs = inspect.getargspec(equations.ADERDG.__init__).args[1:]
cmdArgsDict = vars(cmdLineArgs)
cmdArgsDict['memLayout'] = mem_layout
args = [cmdArgsDict[key] for key in adgArgs]
adg = equations.ADERDG(*args)


compiler = Generator(arch)

# Equation-specific kernels
adg.addInit(compiler)
adg.addLocal(compiler, platforms)
adg.addNeighbor(compiler, platforms)
adg.addTime(compiler, platforms)

# Common kernels
DynamicRupture.addKernels(compiler,
                          adg,
                          cmdLineArgs.matricesDir,
                          cmdLineArgs.dynamicRuptureMethod,
                          platforms)

Plasticity.addKernels(compiler, adg, cmdLineArgs.matricesDir, cmdLineArgs.PlasticityMethod)
SurfaceDisplacement.addKernels(compiler, adg)
Point.addKernels(compiler, adg)

gemm_tool_list = cmdLineArgs.gemm_tools.replace(" ", "").split(",")
generators = []
for tool in gemm_tool_list:
  if hasattr(gemm_configuration, tool):
    specific_gemm_class = getattr(gemm_configuration, tool)
    generators.append(specific_gemm_class(arch))
  else:
    raise RuntimeError(f'YATETO::ERROR: unknown \"{tool}\" GEMM tool. '
                       f'Please, refer to the documentation')

gemmTools = GeneratorCollection(generators)

# Generate code
compiler.generate(cmdLineArgs.outputDir, 'seissol', gemmTools)
