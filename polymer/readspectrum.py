with open("water.out", "r") as f:
    lines = f.readlines()
    result = False
    result_ls = []
    for idx, line in enumerate(lines, 1):
        if "Results With intensities" in line:
            result = not result
        if result:
            line = line.strip().split(" ")
            line = [elem for elem in line if elem.strip()]
            result_ls.append(line)

for i in result_ls[5:-3]:
    print(i[1], i[3])

