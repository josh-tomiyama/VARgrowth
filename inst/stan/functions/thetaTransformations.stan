// // This file contains convenient function transformations for parameters
//
real inverseLogit(real x){
  return(exp(x)/(1 + exp(x)));
}

vector inverseLogitVec(vector x){
  return(exp(x) ./ (rep_vector(1, num_elements(x)) + exp(x)) );
}



vector transformThetaVec(vector transformType, vector theta){
  vector[num_elements(theta)] transform_theta;

  if(num_elements(transformType) != num_elements(theta)){
    reject("Number of transformations not equal to number of theta parameters");
  }

  for(i in 1:num_elements(theta)){
    if(transformType[i] == 1){
      transform_theta[i] = theta[i];
    }else if(transformType[i] == 2){
      transform_theta[i] = exp(theta[i]);
    }else if(transformType[i] == 3){
      transform_theta[i] = inverseLogit(theta[i]);
    }else{
      reject("transform type is not valid");
    }
  }
  return(transform_theta);
}


matrix transformThetaMat(vector transformType, matrix theta){
  matrix[rows(theta), cols(theta)] transform_theta;

  if(transformType[1] == 0){
    transform_theta = UserDefinedThetaTransformation(theta);
    return(transform_theta);
  }

  if(num_elements(transformType) != cols(theta)){
    reject("Number of transformations not equal to number of theta parameters");
  }

  for(i in 1:cols(theta)){
    if(transformType[i] == 1){
      transform_theta[,i] = theta[,i];
    }else if(transformType[i] == 2){
      transform_theta[,i] = exp(theta[,i]);
    }else if(transformType[i] == 3){
      transform_theta[,i] = inverseLogitVec(theta[,i]);
    }else{
      reject("transform type is not valid");
    }
  }
  return(transform_theta);
}

 matrix expLogitLinkTheta(matrix theta){
  matrix[rows(theta), cols(theta)] transform_theta;
  transform_theta[:,1] = exp(theta[:,1]);
  transform_theta[:,2] = exp(theta[:,2]);
  transform_theta[:,3] = exp(theta[:,3]) ./ (rep_vector(1, rows(theta)) + exp(theta[:,3]));
  return(transform_theta);
}
