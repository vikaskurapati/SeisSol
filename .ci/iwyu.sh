iwyu_tool.py -p compile_commands.json -j $(nproc) ../src/ > iwyu.out
fix_includes.py -b --nocomments -m < iwyu.out
