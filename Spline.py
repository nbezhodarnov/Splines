import numpy as np
import matplotlib.pyplot as plt
import math
	
class Spline():
	x = np.array([], dtype = float)
	y = np.array([], dtype = float)
	n = 0
	y_ = None
	coefficients = None
	spline_type = 0
	
	def __init__(self, x_array_input, y_array_input, count, type_of_spline, derivatives):
		self.x = x_array_input
		self.y = y_array_input
		self.n = count
		self.spline_type = type_of_spline
		if (type_of_spline == 1):
			self.coefficients = np.zeros((2, count - 1), dtype = float)
		elif (type_of_spline == 2):
			self.coefficients = np.zeros((4, count - 1), dtype = float)
		elif (type_of_spline == 3):
			self.coefficients = np.zeros((4, count - 1), dtype = float)
			self.y_ = derivatives
		self.__Coefficients_calculate()
	
	def __Coefficients_calculate(self):
		if (self.spline_type == 1):
			for i in range(self.n - 1):
				self.coefficients[0][i] = self.y[i]
				self.coefficients[1][i] = (self.y[i + 1] - self.y[i]) / (self.x[i + 1] - self.x[i])
			return
		elif (self.spline_type == 2):
			self.coefficients[2][0] = 0
			neww_array = np.zeros((self.n - 1, 1))
			neww_array[0] = 0
			ksi_array = np.zeros((self.n - 1, 1))
			ksi_array[0] = 0
			for i in range(1, self.n - 1):
				ksi_array[i] = (self.x[i] - self.x[i + 1]) / ((self.x[i] - self.x[i - 1]) * ksi_array[i - 1] + 2 * (self.x[i + 1] - self.x[i - 1]))
				neww_array[i] = 3 * ((self.y[i + 1] - self.y[i]) / (self.x[i + 1] - self.x[i]) - (self.y[i] - self.y[i - 1]) / (self.x[i] - self.x[i - 1]) - (self.x[i] - self.x[i - 1]) * neww_array[i - 1]) / ((self.x[i] - self.x[i - 1]) * ksi_array[i - 1] + 2 * (self.x[i + 1] - self.x[i - 1]))
			self.coefficients[2][self.n - 2] = neww_array[self.n - 2]
			for i in range(self.n - 3):
				self.coefficients[2][self.n - i - 3] = ksi_array[self.n - i - 3] * self.coefficients[2][self.n - i - 2] + neww_array[self.n - i - 3]
			for i in range(self.n - 2):
				self.coefficients[0][i] = self.y[i]
				self.coefficients[1][i] = (self.y[i + 1] - self.y[i]) / (self.x[i + 1] - self.x[i]) - (self.x[i + 1] - self.x[i]) * (self.coefficients[2][i + 1] + 2 * self.coefficients[2][i]) / 3
				self.coefficients[3][i] = (self.coefficients[2][i + 1] - self.coefficients[2][i]) / (3 * (self.x[i + 1] - self.x[i]))
			self.coefficients[0][self.n - 2] = self.y[self.n - 2]
			self.coefficients[1][self.n - 2] = (self.y[self.n - 1] - self.y[self.n - 2]) / (self.x[self.n - 1] - self.x[self.n - 2]) - (self.x[self.n - 1] - self.x[self.n - 2]) * 2 * self.coefficients[2][self.n - 2] / 3
			self.coefficients[3][self.n - 2] = -self.coefficients[2][self.n - 2] / (3 * (self.x[self.n - 1] - self.x[self.n - 2]))
		elif (self.spline_type == 3):
			for i in range(self.n - 1):
				self.coefficients[0][i] = self.y[i]
				self.coefficients[1][i] = self.y_[i]
				self.coefficients[2][i] = (3 * self.y[i + 1] - 3 * self.y[i] - 2 * (self.x[i + 1] - self.x[i]) * self.y_[i] - (self.x[i + 1] - self.x[i]) * self.y_[i + 1]) / ((self.x[i + 1] - self.x[i]) ** 2)
				self.coefficients[3][i] = (2 * self.y[i] - 2 * self.y[i + 1] + (self.x[i + 1] - self.x[i]) * self.y_[i] + (self.x[i + 1] - self.x[i]) * self.y_[i + 1]) / ((self.x[i + 1] - self.x[i]) ** 3)
		
	def Spline_calculate(self, x_input):
		index = 1
		while (x_input > self.x[index]):
			index += 1
		if (x_input == self.x[index]):
			return self.y[index]
		index -= 1
		if (self.spline_type == 1):
			return self.coefficients[0][index] + self.coefficients[1][index] * (x_input - self.x[index])
		elif (self.spline_type == 2):
			return self.coefficients[0][index] + self.coefficients[1][index] * (x_input - self.x[index]) + self.coefficients[2][index] * ((x_input - self.x[index]) ** 2) + self.coefficients[3][index] * ((x_input - self.x[index]) ** 3)
		elif (self.spline_type == 3):
			return self.coefficients[0][index] + self.coefficients[1][index] * (x_input - self.x[index]) + self.coefficients[2][index] * ((x_input - self.x[index]) ** 2) + self.coefficients[3][index] * ((x_input - self.x[index]) ** 3)
		return 0
	
def derivative(x_input):
	return -math.sin(x_input) - (2 ** (0.1 * x_input)) * math.log(2) * 0.1
	
def main():
	x = np.array([0, 1.75, 3.5, 5.25, 7])
	y = np.array([0, -1.307, -2.211, -0.927, -0.871])
	y_ = np.array([derivative(0), derivative(1.75), derivative(3.5), derivative(5.25), derivative(7)])
	x_dot = 2.555
	polynomial = Spline(x, y, x.size, 1, y_)
	dots_count = 100
	x_plot_table = np.linspace(0, 7, dots_count, dtype = float)
	y_plot_table = np.linspace(0, 0, dots_count, dtype = float)
	y_plot_table_of_original = np.linspace(0, 0, dots_count, dtype = float)
	for i in range(dots_count):
		y_plot_table[i] = polynomial.Spline_calculate(x_plot_table[i])
		y_plot_table_of_original[i] = math.cos(x_plot_table[i]) - 2 ** (0.1 * x_plot_table[i])
	plt.subplot(221)
	plt.plot(x_plot_table, y_plot_table, 'b-', label = 'Linear Spline')
	plt.plot(x_plot_table, y_plot_table_of_original, 'm--', label = 'Original')
	plt.plot(x, y, 'r*', x_dot, polynomial.Spline_calculate(x_dot), 'o')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.annotate('f(x) ~ ' + str(polynomial.Spline_calculate(x_dot)), xy=(x_dot, polynomial.Spline_calculate(x_dot)), xytext=(x_dot, polynomial.Spline_calculate(x_dot) - 0.32),
             arrowprops=dict(facecolor='black', shrink=0.05),
             )
	plt.legend(loc='upper right')
	
	polynomial = Spline(x, y, x.size, 2, y_)
	y_plot_table = np.linspace(0, 0, dots_count, dtype = float)
	for i in range(dots_count):
		y_plot_table[i] = polynomial.Spline_calculate(x_plot_table[i])
	plt.subplot(222)
	plt.plot(x_plot_table, y_plot_table, 'b-', label = 'Cubic Spline')
	plt.plot(x_plot_table, y_plot_table_of_original, 'm--', label = 'Original')
	plt.plot(x, y, 'r*', x_dot, polynomial.Spline_calculate(x_dot), 'o')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.annotate('f(x) ~ ' + str(polynomial.Spline_calculate(x_dot)), xy=(x_dot, polynomial.Spline_calculate(x_dot)), xytext=(x_dot, polynomial.Spline_calculate(x_dot) - 0.32),
             arrowprops=dict(facecolor='black', shrink=0.05),
             )
	plt.legend(loc='upper right')
	
	polynomial = Spline(x, y, x.size, 3, y_)
	x_plot_table = np.linspace(0, 7, dots_count, dtype = float)
	for i in range(dots_count):
		y_plot_table[i] = polynomial.Spline_calculate(x_plot_table[i])
	plt.subplot(223)
	plt.plot(x_plot_table, y_plot_table, 'b-', label = 'Cubic Hermite Spline')
	plt.plot(x_plot_table, y_plot_table_of_original, 'm--', label = 'Original')
	plt.plot(x, y, 'r*', x_dot, polynomial.Spline_calculate(x_dot), 'o')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.annotate('f(x) ~ ' + str(polynomial.Spline_calculate(x_dot)), xy=(x_dot, polynomial.Spline_calculate(x_dot)), xytext=(x_dot, polynomial.Spline_calculate(x_dot) - 0.32),
             arrowprops=dict(facecolor='black', shrink=0.05),
             )
	plt.legend(loc='upper right')
	
	plt.show()
	
if __name__ == '__main__':
    main()
