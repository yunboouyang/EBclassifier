# EBclassifier
Simulation study code and results for author's manusript *A Empirical Bayes Approach for High Dimensional Classification*.

We proposed 2 linear classifiers: DP classifier and Sparse DP classifier based on nonparametric Bayes estimator proposed in 
*A Nonparametric Bayes Approach for Sparse Sequence Estimation*. Another classifier called Hard Thresh DP classifier is included just for 
comparison. 


We compared our methods with the following well-established methods:
- Empirical Bayes Classifier proposed by Greenshtein and Park
- Independence Rule proposed by Bickel and Levina
- Feature Annealed Independence Rule proposed by Fan and Fan
- Logistic regression adding Lasso penalty

We conducted three simulation studies. The first two simulation studies have the same setting as Greenshtein and Park (2009).
The third one has the same setting as Fan and Fan (2008). In RData file misclassification rates are summarized. We have set seed
for all the simulation studies. A real data example code is also included.
