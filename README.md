# Fractional-Queues
The files presented are codes used in analysis for the paper "Queuing models with Mittag-Leffler inter-event times" (See the link https://doi.org/10.1007/s13540-023-00161-4). The code deals with different models for a single server queue with first-in-first-out discipline where the service/arrival times are distributed according to Mittag-Leffler distributions, which are modelled using the "MittagLeffleR" package in R (See the link https://strakaps.github.io/MittagLeffleR/index.html).

The File "GG1.R" is used in Model 3 of the paper, where the service times and arrival times are given by separate sequences of i.i.d Mittag-Leffler random variables. The code provides a direct simulation of the queue, and at the end the queue length is plotted as a function of time. 

The file "InverseSub.R" is used to plot the trajectory of the inverse stable subordinator with respect to time. These inverse subordinators are used to define Model 1 of the paper, they are also used in the expression of the limit theorem for Model 3. The first two parts of the code construcks two stable subordinators using the "Chambers-Mallows-Stuck" method for simulating stable random variables (see the link https://doi.org/10.1080/01621459.1976.10480344) and then inverts these to construct two inverse subordinators. The last part models the last expression in Theorem 3.3 of our paper, the case where the scaling parameters of arrivals and services match and the limit is expressed as the difference of two inverse subordinators reflected at the origin.

The file "RenewalProb.R" 
