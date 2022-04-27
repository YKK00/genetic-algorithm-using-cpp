# genetic-algorithm-using-cpp
*This folder shows how I use genenic algorithm to solve specific problemes.

Genetic algorithm is a global optimization algorithm, which can find solutions suitable for practical problems by simulating the law of "survival of the fittest" in 
nature. Premature maturity that leads to local convergence and other problems are the main shortcomings of traditional genetic algorithms. This paper adopts the retained optimal descendant algorithm, catastrophe algorithm, stretch fitness algorithm, and adaptive algorithm to improve the traditional genetic algorithm and combine them to
apply to the pre-stack AVO nonlinear inversion problem. 

The inversion application experiments of the classic geological models in China and Mexico have been carried out. The experimental results show that the inversion errors of sensitive parameters such as the rate of change of longitudinal wave velocity and the rate of density change are all controlled within 10%. The improved genetic algorithm in the pre-stack AVO has good applicability to the geological model, which verifies the application value of the algorithm.

By inputting the reflection results in the program, we can use the genetic algorithm to obtain the accurate input value. Compared with the traditional algorithm, the genetic algorithm has the characteristics of prematurely converging to the local extreme value. Therefore, this algorithm successively combines the simulated annealing algorithm, automatic Adaptation algorithm and other algorithms to improve the accuracy of the algorithm.

The error between the result calculated by this algorithm and the real value (various geological parameters) is shown in the following figure:

![GAresults](https://user-images.githubusercontent.com/99771239/165602049-f7d9ec3d-6886-4a32-b02b-a2f79ac6ee3b.jpg)

It can be seen from the inversion results that this algorithm can generate inversion results with high accuracy.
