// User defined functions for growth mean/parameter transformations

vector UserDefinedGrowth(int n_data, int[] groups, real[] times, matrix theta){
  vector[n_data] out;

  for(idx in 1:n_data){
    out[idx] = theta[groups[idx], 1] *
    exp(-theta[groups[idx], 2] *
      theta[groups[idx], 3] ^ times[idx]
    );
  }
  return(out);
}
