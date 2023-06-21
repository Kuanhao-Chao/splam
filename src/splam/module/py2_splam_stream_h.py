f = open("splam_stream.py", "r")

fw = open("splam_stream.h", "w")

with open("splam_stream.py", "r") as f:
    lines = f.read().splitlines()
    fw.write("#include <string>\n")
    fw.write("std::string python_script = \n")
    counter = 0
    for line in lines:
        # print(len(line), line)
        if len(line) > 1:
            print(line)
            if counter < len(lines)-1:
                fw.write('"' + line + '\\n"\n')
            else:
                fw.write('"' + line + '\\n";\n')
        counter += 1
