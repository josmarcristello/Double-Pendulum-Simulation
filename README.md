# A simulation of the double pendulum chaotic motion using Python.

## Introduction
The simple pendulum has an apparent simple motion, and yet, the union of two pendulums generate a system with a chaotic motion. Although there is no formal proof, there is no known solution for the double pendulum system motion, and thus, it must be solved numerically. This study aims to analyze such a system and compare the sensitivity to different numerical methods. Moreover, it is also a goal to goal to analyze how the chaotic motion takes place, and whether it is prevalent for the entire system or for some cases.
To achieve that, equations of motion are obtained using the Euler-Lagrangian equation. The resulting differential equations are them solved numerically, using Runge-Kutta (4th order) and Euler’s Explicit Forward methods. The system is then computationally simulated using the programming language python. Detailed solutions are initially presented for different initial states of the system and the two different numerical methods, analyzing motion, behavior, and mechanical energy of the system. Then, solutions are presented for a very large number of initial positions, allowing for patterns to be discovered.

The notebook in this repository is contains the modelling of the equations of motion, and numerical solution (Implementaiton of the numerical methods). All models used and references can be found on the report and on the [Acknowledgements](#Acknowledgements) and [Bibliography](#Bibliography) sections. The repository is fully open-sourced under the [MIT License](https://choosealicense.com/licenses/mit/). This study was originally submitted as a project for the class ENME631 (Numerical Methods) in the [University of Calgary](https://www.ucalgary.ca/) (Calgary, Alberta, Canada).

![Alt text](./2-figures/1. Numerical Methods Verification/Number of Flips x Initial Angular Positions.png?raw=true "Title")
![Alt text](./2-figures/1. Numerical Methods Verification/Number of Flips x Initial Angular Positions.png)


## Table of Contents
1. [Installation](#Installation)
2. [Model](#Model)
3. [Usage](#Usage)
4. [File Descriptions](#File-Descriptions)
5. [Results](#Results)
6. [Acknowledgements](#Acknowledgements)
7. [Bibliography](#Bibliography)


## Installation
All necessary libraries are shown in the beginning of each notebook. The code was tested with [Python 3.9.7 64-bit](https://www.python.org/downloads/release/python-397/), but should work with version 3.x. 

The recommended way to install the libraries is with [Pypi](https://pypi.org), using pip:
```bash
> pip install <libraryname>
```
Before opening an issue, make sure your libraries are updated:
```bash
> pip install --upgrade <libraryname>
```

## Usage
Either: 

* Clone the repository, which will give you all the data used in the original study. Please see [Acknowledgements](#Acknowledgements) for licensing of each data source information.

* Use the notebooks, which will download and process all the data, and is throughout commented to allow for some flexibility in the simulation.

* It should be noted that the execution time rises exponentially for a large number of points. For example, simulating 1500 initial angular positions for each pendulum bob generates a total of 2,250,000 trajectories, which took approximately 16h of continuous running in the author's machine. A suggestion to improve that is to either parallelize the code or use an implementation such as [Cython](https://cython.org/)

## File Descriptions

1.  [**Double_Pendulum_Model.ipynb**](./1-source/Double_Pendulum_Model.ipynb) - The python notebook with all the models, simulation and numerical solution.
2.  [**Figures Folder**](./2-figures/) - Contains all the images and animations used in the final report. 
2.1 [**Numerical Methods Verification**](./2-figures/1. Numerical Methods Verification/)
2.1 [**Cartesian Motion**](./2-figures/2-Cartesian Motion/)
3.  [**TCC_PUCMG_Report.pdf**](./3-report/ENME 631 - Double Pendulum Motion Analysis.pdf)

## Results
The selected model used Gradient Boosting (XGBoost Library) and had a R² of 0.691 with MSE of 0.125. The notebook had a very extensive data analysis, and contains recommendations for a host that wants to maximize their profit, or a guest that wants to maximize their cost-benefit when renting a place. For the guest example, in Toronto:

## Acknowledgements

* [NumPy.com](https://numpy.org/), the library used for matrix and arrays operation in this project. This library is open-sourced under the [BSD License](https://choosealicense.com/licenses/0bsd/).

* [Matplotlib](https://matplotlib.org/), the library for data visualization. This library uses [BSD compatible code](https://choosealicense.com/licenses/0bsd/).
 and operates in the LSF [License] (https://docs.python.org/3/license.html).

## Bibliography

1. Acheson, D. J. (1997). From calculus to chaos: An introduction to dynamics. Oxford University Press.
2. Stachowiak, T., & Okada, T. (2006). A numerical analysis of chaos in the double pendulum. Chaos, Solitons & Fractals, 29(2), 417–422. https://doi.org/10.1016/j.chaos.2005.08.032
3. Ohlhoff, A., & Richter, P. H. (2000). Forces in the Double Pendulum. ZAMM, 80(8), 517–534. https://doi.org/10.1002/1521-4001(200008)80:8<517::AID-ZAMM517>3.0.CO;2-1
4. Marcelo Tusset, A., Piccirillo, V., Bueno, A. M., Manoel Balthazar, J., Sado, D., Felix, J. L. P., & Brasil, R. M. L. R. da F. (2016). Chaos control and sensitivity analysis of a double pendulum arm excited by an RLC circuit based nonlinear shaker. Journal of Vibration and Control, 22(17), 3621–3637. https://doi.org/10.1177/1077546314564782
5. Sanjeewa, S. D., & Parnichkun, M. (2019). Control of rotary double inverted pendulum system using mixed sensitivity H∞ controller. International Journal of Advanced Robotic Systems, 16(2), 172988141983327. https://doi.org/10.1177/1729881419833273
6. Kiyoumarsi, A., Ataei, M., Mirzaeian-Dehkordi, B., & Ashrafi, R. (2007). The Mathematical Modeling of a Double-Pendulum System as a Physical Model of Flexible Arm Robot. 2007 IEEE International Conference on Control and Automation, 1900–1904. https://doi.org/10.1109/ICCA.2007.4376692
7. Sun, N., Wu, Y., Liang, X., & Fang, Y. (2019). Nonlinear Stable Transportation Control for Double-Pendulum Shipboard Cranes With Ship-Motion-Induced Disturbances. IEEE Transactions on Industrial Electronics, 66(12), 9467–9479. https://doi.org/10.1109/TIE.2019.2893855
8. Bazargan-Lari, Y., Eghtesad, M., Khoogar, A., & Mohammad-Zadeh, A. (2014). Dynamics and regulation of locomotion of a human swing leg as a double-pendulum considering self-impact joint constraint. Journal of Biomedical Physics & Engineering, 4(3), 91–102.
9. Inoue, Y., Shibata, K., Fukumoto, M., Wang, S., & Oka, K. (n.d.). Study on Dynamic Analysis and Wearable Sensor System for Golf Swing. 146.
10. Neal, H. (n.d.). Using Mechanics of a Double Pendulum to Maximize Sport Performance. 6.
11. Gonzalez, G. (2008). Single and Double plane pendulum. 8.
12. Zingg, D. W., & Chisholm, T. T. (1999). Runge–Kutta methods for linear ordinary differential equations. Applied Numerical Mathematics, 31(2), 227–238. https://doi.org/10.1016/S0168-9274(98)00129-9
13. J, R. (1983). Runge-kutta method. Numerical Methods.
14. Levien, R. B., & Tan, S. M. (1998). Double pendulum: An experiment in chaos. American Journal of Physics, 61(11), 1038. https://doi.org/10.1119/1.17335
15. Chen, Joe. (2008). Chaos from Simplicity: An Introduction to the Double Pendulum.

 

