// User defined functions for growth mean/parameter transformations

matrix UserDefinedThetaTransformation(matrix theta){
  matrix[rows(theta), cols(theta)] out;
  out = theta;
  return(out);
}
