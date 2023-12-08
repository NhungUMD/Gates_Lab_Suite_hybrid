# Daiwei Zhu Sep-26-2018
# Basic definition of Quantum Gates/Quantum States/Quantum Circuits

# Please Do not circulate

import random
import os
# from scipy.stats import binom
from scipy.stats import *
import scipy.linalg as linalg
import numpy as np
__author__ = "Daiwei Zhu"


global param_table
param_table = dict()

global s
s = [[[1, 0], [0, 1]], [[0, 1], [1, 0]], [
    [0, -1J], [1J, 0]], [[1, 0], [0, -1]]]

global mapping
mapping = [i for i in range(1, 33)]

global mapping_back
mapping_back = [i for i in range(-1, 32)]

global pi
pi = 3.141592653589793

global max_size
max_size = 7


def Identity(size):
    matrix = 1
    for i in range(1, size+1):
        matrix = np.kron(matrix, s[0])
    return matrix


def XX_Rotation(size, qubit1, qubit2, param):
    matrix = -1J*param
    for i in range(size):
        if (qubit1 == i) or (qubit2 == i):
            matrix = np.kron(matrix, s[1])
        else:
            matrix = np.kron(matrix, s[0])
    return linalg.expm(matrix)


def X_Rotation(size, qubit, param):
    matrix = -1J*param/2
    for i in range(size):
        if qubit == i:
            matrix = np.kron(matrix, s[1])
        else:
            matrix = np.kron(matrix, s[0])
    return linalg.expm(matrix)


def Y_Rotation(size, qubit, param):
    matrix = -1J*param/2
    for i in range(size):
        if qubit == i:
            matrix = np.kron(matrix, s[2])
        else:
            matrix = np.kron(matrix, s[0])
    return linalg.expm(matrix)


def Z_Rotation(size, qubit, param):
    matrix = -1J*param/2
    for i in range(size):
        if qubit == i:
            matrix = np.kron(matrix, s[3])
        else:
            matrix = np.kron(matrix, s[0])
    return linalg.expm(matrix)


# def inver_bino_lower(x, n, x1, x2):
#     temp = (x1+x2)/2
#     if abs(temp-x1) < 0.000001:
#         return temp
#     elif binom.cdf(x, n, temp) > 0.84135:
#         return inver_bino_upper(x, n, temp, x2)
#     else:
#         return inver_bino_upper(x, n, x1, temp)


# def inver_bino_upper(x, n, x1, x2):
#     """
#     Determines the value 1 standard deviation above the mean for a binomial distribution

#     params:
#         x: number of bright counts (successes)
#         n: number of trials
#         x1: ?
#         x2: ?
#     """
#     temp = (x1+x2)/2
#     if abs(temp-x1) < 0.000001:
#         return temp
#     elif binom.cdf(x, n, temp) > 0.15865:
#         return inver_bino_lower(x, n, temp, x2)
#     else:
#         return inver_bino_lower(x, n, x1, temp)


# def error(x, n):
#     """
#     Determines the upper and lower values 1 standard deviation from the mean for a binomial distribution

#     # I DON'T THINK THIS WORKS - CHECK IT WITH MATLAB VERSION

#     parameters:
#         x: bright counts
#         n: total trials
#     """
#     return inver_bino_lower(x, n, 0, 1), inver_bino_upper(x, n, 0, 1)


def error(x, n, alpha=0.3173):
    """
    Produces lower and upper bounds from a binomial distribution.

    params:
        x: number of successes
        n: number of total trials
        alpha: confidence interval (default of 0.15865 * 2 for one sigma)
    """

    # Lower limit:
    nu1 = 2 * x
    nu2 = 2 * (n - x + 1)
    F = f.ppf(alpha / 2, nu1, nu2)
    lb = (nu1 * F) / (nu2 + nu1 * F)

    # Upper limit:
    nu1 = 2 * (x + 1)
    nu2 = 2 * (n - x)
    F = f.ppf(1 - alpha / 2, nu1, nu2)
    ub = (nu1 * F) / (nu2 + nu1 * F)

    return lb, ub


def Define_Parameter(name, value=0):
    param_table[name] = value


def Get_Parameter(sequence):
    temp = sequence.split("*")
    value = 1
    for key in temp:
        value = value*param_table[key]
    return value


def Regulate_Params():
    for key in param_table:
        param_table[key] = param_table[key] % 2
        if param_table[key] >= 1:
            param_table[key] -= 2
        if param_table[key] <= -1:
            param_table[key] += 2


def Energy_Z(state):
    temp_state = np.multiply(np.add(state, -0.5), -2)
    return -np.sum(temp_state)


def Correlator(i, j):
    return lambda state: (state[i]-0.5)*(state[j]-0.5)*4


def Qubit(i):
    return lambda state: (state[i]-0.5)*(-2)


def Energy_ZZ(state):
    temp_state = np.multiply(np.add(state, -0.5), -2)
    counter = 0
    for i in range(len(temp_state)-1):
        counter += temp_state[i]*temp_state[i+1]
    counter += temp_state[0]*temp_state[len(temp_state)-1]
    return -counter


def Binary_Encoding(data_array):
    code = 0
    for i in range(len(data_array)):
        code = code*2+data_array[i]
    return code


def Binary_Decoding(code, digits):
    data_array = [0]*digits
    counter = code
    for i in reversed(range(digits)):
        data_array[i] = counter % 2
        counter = counter//2
    return data_array


def Reset_Mapping():
    mapping = range(1, 33)
    mapping_back = range(-1, 32)


def Set_Mapping(index):
    for i in range(len(index)):
        mapping[i] = index[i]
    for i in range(len(index)):
        mapping_back[index[i]] = i


def Trace_Out(state, indicator):
    temp = Quantum_State(len(indicator))
    temp_index = [0]*temp.size
    for i in range(len(temp.population)):
        temp.population[i] = 0

    for i in range(len(state.population)):
        index = Binary_Decoding(i, state.size)
        for j in range(temp.size):
            temp_index[j] = index[indicator[j]]
        code = Binary_Encoding(temp_index)
        temp.population[code] += state.population[i]

    assert (abs(np.sum(temp.population)-1) <
            0.0000001), "population not normalized"

    return temp


class Quantum_State:
    detect_threshold = 1

    def __init__(self, size):
        self.size = size
        self.state = [0]*(2**size)
        self.state[0] = 1
        self.population_only = False
        self.population = [0]*(2**size)
        self.error_lower = [0]*(2**size)
        self.error_upper = [0]*(2**size)
        self.Update_Population()

    def Normalization(self):
        temp_sum = np.sum(self.population)
        self.population = np.divide(self.population, temp_sum)
        if self.population_only == True:
            return
        temp_sum = np.sqrt(temp_sum)
        self.state = np.divide(self.state, temp_sum)
        return

    def ImportFromList(self, data_list):
        """
        Get population data from raw data list
        """
        self.population_only = True
        self.population = np.multiply(self.population, 0)

        for line in data_list:
            temp = list(map(int, line))[1:]
            current_line = [0]*self.size
            for i in range(self.size):
                if temp[mapping[i]] > self.detect_threshold:
                    current_line[i] = 1
                else:
                    current_line[i] = 0

            code = Binary_Encoding(current_line)
            self.population[code] += 1

        temp_sum = np.sum(self.population)

        for i in range(len(self.population)):
            self.error_lower[i], self.error_upper[i] = error(
                self.population[i], temp_sum)

        self.population = np.divide(self.population, temp_sum)

    def Import(self, filename, **kwarg):
        """
        Get Population data from raw data files
        """
        read_file = open(filename, "r")
        self.population_only = True

        # reset populations all to 0?
        self.population = np.multiply(self.population, 0)

        # Not great coding practice - should probably be changed eventually
        # Go through each line of the file
        while True:
            temp = read_file.readline().split()
            if len(temp) == 0:
                break
            # what is this
            temp = list(map(int, temp))[1:]
            current_line = [0]*self.size
            for i in range(self.size):
                if temp[mapping[i]] > self.detect_threshold:
                    current_line[i] = 1
                else:
                    current_line[i] = 0
            code = Binary_Encoding(current_line)
            self.population[code] += 1
        temp_sum = np.sum(self.population)
        for i in range(len(self.population)):
            self.error_lower[i], self.error_upper[i] = error(
                self.population[i], temp_sum)
        self.population = np.divide(self.population, temp_sum)
        read_file.close()

        if "keep_file" in kwarg:
            if kwarg["keep_file"] == True:
                return

        os.remove(filename)

    def Mean_Value(self, function):
        mean_value = 0
        for i in range(len(self.population)):
            mean_value += self.population[i] * \
                function(Binary_Decoding(i, self.size))
        return mean_value

    def Std(self, function):
        std = 0
        for i in range(len(self.population)):
            std += np.square((self.error_upper[i]-self.error_lower[i])/2)*np.square(
                function(Binary_Decoding(i, self.size)))
        return np.sqrt(std)

    def Update_Population(self):
        if self.population_only:
            raise ErrorValue("state doesn't contain phase information")
        for i in range(len(self.state)):
            self.population[i] = abs(self.state[i])**2

    def Measure(self, file_name, sampling_number=2000):
        random.seed()
        write_file = open(file_name, 'w')
        for i in range(sampling_number):
            j = random.random()
            for k in range(len(self.population)):
                j -= self.population[k]
                if j <= 0:
                    break
            qubit_state = Binary_Decoding(k, self.size)
            ion_state = [0]*(max_size+1)
            for k in range(self.size):
                ion_state[mapping[k]] = qubit_state[k]*5
            to_write_string = "0\t0"
            for k in range(1, max_size+1):
                to_write_string = to_write_string+"\t"+str(int(ion_state[k]))
            write_file.write(to_write_string+"\n")
        write_file.close()

        return


class Quantum_Gate:
    def __init__(self, name, qubit1, qubit2=None, **kwarg):
        self.name = name
        self.qubit1 = qubit1
        self.qubit2 = qubit2
        if "angle" in kwarg:
            self.angle = kwarg["angle"]
        else:
            self.angle = None

        if "axis" in kwarg:
            self.axis = kwarg["axis"]
        else:
            self.axis = None

    def Matrix_Representation(self, size):
        if self.angle != None:
            try:
                param = float(self.angle)
            except:
                param = param_table[self.angle]

        if self.axis != None:
            try:
                axis = float(self.axis)
            except:
                axis = param_table[self.axis]

        if self.name == "FTXA":
            return XX_Rotation(size, self.qubit1, self.qubit2, param*pi)

        elif self.name == "XA":
            return XX_Rotation(size, self.qubit1, self.qubit2, param*pi)

        elif (self.name == "AZ"):
            return Z_Rotation(size, self.qubit1, param*pi)

        elif (self.name == "SK1Z"):
            return Z_Rotation(size, self.qubit1, pi/param)

        elif (self.name == "SKAX"):
            return X_Rotation(size, self.qubit1, param*pi)

        elif (self.name == "SK1X"):
            return X_Rotation(size, self.qubit1, pi/param)

        elif (self.name == "SKAY"):
            return Y_Rotation(size, self.qubit1, pi*param)

        elif (self.name == "SKAAZY"):
            return Y_Rotation(size, self.qubit1, pi*param)

        elif (self.name == "SK1Y"):
            return Y_Rotation(size, self.qubit1, pi/param)

        elif (self.name == "FTXX"):
            if (self.angle != None):
                return XX_Rotation(size, self.qubit1, self.qubit2, pi*param)
            else:
                return XX_Rotation(size, self.qubit1, self.qubit2, pi/4)

        elif (self.name == "XX"):
            if self.angle != None:
                return XX_Rotation(size, self.qubit1, self.qubit2, pi*param)
            else:
                return XX_Rotation(size, self.qubit1, self.qubit2, pi/4)

        elif (self.name == "CYA"):
            temp = Z_Rotation(size, self.qubit1, -0.5*pi)
            temp = np.matmul(Z_Rotation(size, self.qubit2, -0.5*pi), temp)
            temp = np.matmul(XX_Rotation(size, self.qubit1,
                             self.qubit2, param*pi), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit1, 0.5*pi), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit2, 0.5*pi), temp)
            return temp

        elif (self.name == "FTCYA"):
            temp = Z_Rotation(size, self.qubit1, -0.5*pi)
            temp = np.matmul(Z_Rotation(size, self.qubit2, -0.5*pi), temp)
            temp = np.matmul(XX_Rotation(size, self.qubit1,
                             self.qubit2, param*pi), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit1, 0.5*pi), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit2, 0.5*pi), temp)
            return temp

        elif (self.name == "RZ"):
            return Z_Rotation(size, self.qubit1, pi/param)

        elif (self.name == "RX"):
            return X_Rotation(size, self.qubit1, pi/param)

        elif (self.name == "RY"):
            return Y_Rotation(size, self.qubit1, pi/param)

        elif (self.name == "HSB"):
            temp = XX_Rotation(size, self.qubit1, self.qubit2, param*pi)
            temp = np.matmul(Y_Rotation(size, self.qubit1, pi/2), temp)
            temp = np.matmul(Y_Rotation(size, self.qubit2, pi/2), temp)
            temp = np.matmul(XX_Rotation(size, self.qubit1,
                             self.qubit2, param*pi), temp)
            temp = np.matmul(Y_Rotation(size, self.qubit1, -pi/2), temp)
            temp = np.matmul(Y_Rotation(size, self.qubit2, -pi/2), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit1, -pi/2), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit2, -pi/2), temp)
            temp = np.matmul(XX_Rotation(size, self.qubit1,
                             self.qubit2, param*pi), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit1, pi/2), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit2, pi/2), temp)
            return temp
        elif (self.name == "CNOT"):
            temp = Y_Rotation(size, self.qubit1, pi/2)
            temp = np.matmul(XX_Rotation(
                size, self.qubit1, self.qubit2, pi/4), temp)
            temp = np.matmul(Y_Rotation(size, self.qubit1, -pi/2), temp)
            temp = np.matmul(X_Rotation(size, self.qubit2, -pi/2), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit1, -pi/2), temp)
            return temp
        elif (self.name == "CN"):
            temp = Y_Rotation(size, self.qubit1, pi/2)
            temp = np.matmul(XX_Rotation(
                size, self.qubit1, self.qubit2, pi/4), temp)
            temp = np.matmul(Y_Rotation(size, self.qubit1, -pi/2), temp)
            temp = np.matmul(X_Rotation(size, self.qubit2, -pi/2), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit1, -pi/2), temp)
            return temp

        elif (self.name == "SKH"):
            return np.matmul(Y_Rotation(size, self.qubit1, -pi/2), X_Rotation(size, self.qubit1, pi))
        elif (self.name == "HAD"):
            return np.matmul(Y_Rotation(size, self.qubit1, -pi/2), X_Rotation(size, self.qubit1, pi))

        elif (self.name == "RA"):
            temp = Z_Rotation(size, self.qubit1, -pi*axis)
            temp = np.matmul(X_Rotation(size, self.qubit1, pi*param), temp)
            temp = np.matmul(Z_Rotation(size, self.qubit1, pi*axis), temp)
            return temp

        else:
            raise ValueError("Gate is not defined")

    def GatesLab_Command(self):
        if self.angle != None:
            try:
                param = float(self.angle)
            except:
                param = param_table[self.angle]

        if self.axis != None:
            try:
                axis = float(self.axis)
            except:
                axis = param_table[self.axis]

        if self.name.find("SKA") > -1 and self.name.find("SKAAZY") == -1:
            section = self.name+str(mapping[self.qubit1])
            if param > 0:
                section += "+"
            else:
                section += "-"
            section += str("{:.4f}").format(abs(param))
        elif self.name.find("SKAAZY") > -1:
            section = self.name+str(mapping[self.qubit1])
            if param > 0:
                section += "+"
            else:
                section += "-"
        elif self.name.find("SKH") > -1:
            section = self.name+str(mapping[self.qubit1])
        elif self.name.find("HAD") > -1:
            section = self.name+str(mapping[self.qubit1])

        elif self.name.find("CNOT") > -1:
            section = self.name + \
                str(mapping[self.qubit1])+str(mapping[self.qubit2])

        elif self.name.find("CN") > -1:
            section = self.name + \
                str(mapping[self.qubit1])+str(mapping[self.qubit2])

        elif self.name.find("SK1") > -1:
            section = self.name
            if param > 0:
                section += "+"
            else:
                section += "-"
            section += str(mapping[self.qubit1])+str(int(abs(param)))

        elif self.name.find("FTXA") > -1:
            if param >= 0.5:
                param -= 1
            if param <= -0.5:
                param += 1
            section = self.name + \
                str(mapping[self.qubit1])+str(mapping[self.qubit2])
            section += str("{:.4f}").format(abs(param))
            if param > 0:
                section += "+"
            else:
                section += "-"
        elif self.name.find("XA") > -1:
            if param >= 0.5:
                param -= 1
            if param <= -0.5:
                param += 1
            section = self.name + \
                str(mapping[self.qubit1])+str(mapping[self.qubit2])
            section += str("{:.4f}").format(abs(param))
            if param > 0:
                section += "+"
            else:
                section += "-"
        elif self.name.find("CYA") > -1 and self.name.find("FTCYA") == -1:
            if param >= 0.5:
                param -= 1
            if param <= -0.5:
                param += 1
            section = self.name + \
                str(mapping[self.qubit1])+str(mapping[self.qubit2])
            section += str("{:.4f}").format(abs(param))
            if param > 0:
                section += "+"
            else:
                section += "-"

        elif self.name.find("FTCYA") > -1:
            if param >= 0.5:
                param -= 1
            if param <= -0.5:
                param += 1
            section = self.name + \
                str(mapping[self.qubit1])+str(mapping[self.qubit2])
            section += str("{:.4f}").format(abs(param))
            if param > 0:
                section += "+"
            else:
                section += "-"
        elif self.name.find("AZ") > -1:
            section = self.name+str(mapping[self.qubit1])
            if param > 0:
                section += "+"
            else:
                section += "-"
            section += str("{:.4f}").format(abs(param))
        elif self.name.find("FTXX") > -1:
            if self.angle != None:
                if param < 0:
                    section = self.name + \
                        str(mapping[self.qubit1])+str(mapping[self.qubit2])+"-"
            else:
                section = self.name + \
                    str(mapping[self.qubit1])+str(mapping[self.qubit2])+"+"
        elif self.name.find("XX") > -1:
            if self.angle != None:
                if param < 0:
                    section = self.name + \
                        str(mapping[self.qubit1])+str(mapping[self.qubit2])+"-"
            else:
                section = self.name + \
                    str(mapping[self.qubit1])+str(mapping[self.qubit2])+"+"

        elif self.name.find("HSB") > -1:
            if param >= 0.5:
                param -= 1
            if param <= -0.5:
                param += 1
            if mapping[self.qubit1] < mapping[self.qubit2]:
                section = self.name + \
                    str(mapping[self.qubit1])+str(mapping[self.qubit2])
            else:
                section = self.name + \
                    str(mapping[self.qubit2])+str(mapping[self.qubit1])
            if param > 0:
                section += "+"
            else:
                section += "-"
            section += str("{:.4f}").format(abs(param))

        elif self.name.find("RA") > -1:
            section = self.name+str(mapping[self.qubit1])
            if param > 0:
                section += "+"
            else:
                section += "-"
            section += str("{:.4f}").format(abs(param))+","
            section += str("{:.4f}").format(abs(axis))

        elif (self.name.find("RZ") > -1) or (self.name.find("RX") > -1) or (self.name.find("RY") > -1):
            section = self.name
            if param > 0:
                section += "+"
            else:
                section += "-"
            section += str(mapping[self.qubit1])+str(int(abs(param)))

        return section


class Quantum_Circuit:
    def __init__(self, size, name):
        self.size = size
        self.depth = 0
        self.gates = []
        self.name = name

    def Interpret_GatesLab_Sequence(self, sequence):
        # all the ion indices will be remapped into qubit indices at the end of this block.
        self.depth = 0
        self.gates = []
        temp = sequence.split(":")
        for command in temp:
            if command.find("SKA") > -1 and command.find("SKAAZ") == -1:
                name = command[0:4]
                qubit1 = int(command[4:5])
                param = float(command[5:])
                self.Add_Gate(Quantum_Gate(name, qubit1, angle=param))
            elif command.find("SKAAZY") > -1:
                name = command[0:6]
                qubit1 = int(command[6:7])
                param = float(command[7:13])
#                 param2=float(command[13:19])
                self.Add_Gate(Quantum_Gate(name, qubit1, angle=param))

            elif command.find("SK1") > -1:
                name = command[0:4]
                qubit1 = int(command[5:6])
                param = float(command[6:])
                if command.find("-") > -1:
                    param = param*(-1)
                self.Add_Gate(Quantum_Gate(name, qubit1, angle=param))
            elif command.find("FTXX") > -1:
                name = command[0:4]
                qubit1 = int(command[4:5])
                qubit2 = int(command[5:6])
                param = 0.25
                if command.find("-") > -1:
                    param = param*(-1)
                self.Add_Gate(Quantum_Gate(name, qubit1, qubit2, angle=param))
            elif command.find("XX") > -1:
                name = command[0:2]
                qubit1 = int(command[2:3])
                qubit2 = int(command[3:4])
                param = 0.25
                if command.find("-") > -1:
                    param = param*(-1)
                self.Add_Gate(Quantum_Gate(name, qubit1, qubit2, angle=param))

            elif command.find("FTXA") > -1:
                name = command[0:4]
                qubit1 = int(command[4:5])
                qubit2 = int(command[5:6])
                param = float(command[6:12])
                if command.find("-") > -1:
                    param = param*(-1)
                self.Add_Gate(Quantum_Gate(name, qubit1, qubit2, angle=param))

            elif command.find("XA") > -1:
                name = command[0:2]
                qubit1 = int(command[2:3])
                qubit2 = int(command[3:4])
                param = float(command[4:10])
                if command.find("-") > -1:
                    param = param*(-1)
                self.Add_Gate(Quantum_Gate(name, qubit1, qubit2, angle=param))
            elif command.find("CYA") > -1 and command.find("FTCYA") == -1:
                name = command[0:3]
                qubit1 = int(command[3:4])
                qubit2 = int(command[4:5])
                param = float(command[5:11])
                if command.find("-") > -1:
                    param = param*(-1)
                self.Add_Gate(Quantum_Gate(name, qubit1, qubit2, angle=param))
            elif command.find("FTCYA") > -1:
                name = command[0:5]
                qubit1 = int(command[5:6])
                qubit2 = int(command[6:7])
                param = float(command[7:13])
                if command.find("-") > -1:
                    param = param*(-1)
                self.Add_Gate(Quantum_Gate(name, qubit1, qubit2, angle=param))
            elif command.find("AZ") > -1:
                name = command[0:2]
                qubit1 = int(command[2:3])
                param = float(command[3:])
                self.Add_Gate(Quantum_Gate(name, qubit1, angle=param))

            elif command.find("SKH") > -1:
                name = command[0:3]
                qubit1 = int(command[3:4])
                self.Add_Gate(Quantum_Gate(name, qubit1))
            elif command.find("HAD") > -1:
                name = command[0:3]
                qubit1 = int(command[3:4])
                self.Add_Gate(Quantum_Gate(name, qubit1))

            elif command.find("HSB") > -1:
                name = command[0:3]
                qubit1 = int(command[3:4])
                qubit2 = int(command[4:5])
                param = float(command[5:])
                if command.find("-") > -1:
                    param = param*(-1)
                self.Add_Gate(Quantum_Gate(name, qubit1, qubit2, angle=param))

            elif command.find("RA") > -1:
                name = command[0:2]
                qubit1 = int(command[2:3])
                param = float(command[3:10])
                axis = float(command[11:])
                self.Add_Gate(Quantum_Gate(
                    name, qubit1, angle=param, axis=axis))

            elif command.find("CNOT") > -1:
                name = command[0:4]
                qubit1 = int(command[4:5])
                qubit2 = int(command[5:6])
                self.Add_Gate(Quantum_Gate(name, qubit1, qubit2))

            elif command.find("CN") > -1:
                name = command[0:2]
                qubit1 = int(command[2:3])
                qubit2 = int(command[3:4])
                self.Add_Gate(Quantum_Gate(name, qubit1, qubit2))

            elif command.find("RZ") > -1 or command.find("RY") > -1 or command.find("RX") > -1:
                name = command[0:2]
                qubit1 = int(command[3:4])
                param = int(command[4:])
                if command.find("-") > -1:
                    param = param*(-1)
                self.Add_Gate(Quantum_Gate(name, qubit1, angle=param))

            else:
                raise ValueError("cannot identify"+command)

            # mapping the ion index into qubit index
            self.gates[self.depth -
                       1].qubit1 = mapping_back[self.gates[self.depth-1].qubit1]
            if self.gates[self.depth-1].qubit2 != None:
                self.gates[self.depth -
                           1].qubit2 = mapping_back[self.gates[self.depth-1].qubit2]

        if (self.gates[self.depth-1].qubit1 >= self.size) or (self.gates[self.depth-1].qubit2 != None and self.gates[self.depth-1].qubit2 >= self.size):
            raise ValueError(
                "number of qubit involved in the circuit is more than that defined for the circuit.")
        return

    def Matrix_Representation(self):
        matrix = Identity(self.size)
        for gate in self.gates:
            matrix = np.matmul(gate.Matrix_Representation(self.size), matrix)

        return matrix

    def Simulate(self):
        temp = Quantum_State(self.size)
        temp.state = np.matmul(self.Matrix_Representation(), temp.state)
        temp.Update_Population()
        return temp

    def Emulate(self, file_name, sampling_number=2000):
        temp = self.Simulate()
        temp.Measure(file_name, sampling_number)
        return

    def Add_Gate(self, quantum_gate):
        self.depth += 1
        self.gates.append(quantum_gate)

    def GatesLab_Sequence(self):
        sequence = ""
        for i in range(self.depth):
            sequence += self.gates[i].GatesLab_Command()+":"
        return sequence[:-1]

    def GatesLab_Sequence_XMeasurement(self):
        sequence = ""
        for i in range(self.depth):
            sequence += self.gates[i].GatesLab_Command()+":"
        temp_gate = Quantum_Gate("SK1Y", 0, angle="-2")
        for i in range(self.size):
            temp_gate.qubit1 = i
            sequence += temp_gate.GatesLab_Command()+":"
        return sequence[:-1]
