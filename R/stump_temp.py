# fit(X, y, sample_weight=None, seed or self.seed)

# no scaling

loop on predictors in cols: 
  p
  
  loop on cutpoints: (np.unique?) # omp
    p -> cutpoint (numeric) =: s
        -> averages / majority voting on p < s and p >= s
          -> keep avgs / maj. votes
            -> get training error: (weighted) RMSE / (weighted) error rate
    
return (predictor index, cutpoint, avgs/maj. votes) w/ min (weighted) RMSE / (weighted) error rate
