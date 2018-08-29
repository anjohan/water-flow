def extract_variable(varname):
    with open("log.system", "r") as infile:
        lines = infile.readlines()

    searchstring = "variable " + varname + " equal"
    for line in lines[::-1]:
        if searchstring in line:
            result = eval(line.split()[-1])
            return result

    errorstring = "Variable " + varname + " not found!"
    raise Exception(errorstring)


if __name__ == "__main__":
    print(extract_variable("z_topofcylinder"))
