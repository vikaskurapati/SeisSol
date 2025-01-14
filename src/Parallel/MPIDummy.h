/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Dummy MPI Wrapper
 */

#ifndef MPI_DUMMY_H
#define MPI_DUMMY_H

#include "Parallel/MPIBasic.h"

namespace seissol
{

// MPI_Comm substitute, so we don't have to disambiguate when it's not there
#ifndef USE_MPI
using MPI_Comm = int;
constexpr MPI_Comm MPI_COMM_NULL = 0;
constexpr MPI_Comm MPI_COMM_SELF = 1;
constexpr MPI_Comm MPI_COMM_WORLD = 2;
#endif

/**
 * Dummy class when running without MPI
 */
class MPIDummy : public MPIBasic
{
private:
	MPIDummy()
	{ }

public:
	~MPIDummy()
	{ }

	/**
	 * Does nothing
	 */
	void init(int &argc, char** &argv)
	{
	}

	/**
	 * @return Dummy communicator
	 */
	int comm() const
	{
		return 0;
	}

	/**
	 * Does nothing
	 */
	void barrier(int comm) const
	{
	}

	/**
	 * Does nothing
	 */
	void finalize()
	{
	}

	/** The only instance of the class */
	static MPIDummy mpi;
};

}

#endif // MPI_DUMMY_H
