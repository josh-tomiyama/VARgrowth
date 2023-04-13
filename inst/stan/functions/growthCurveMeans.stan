// //
// // Series of convenience functions to help with VAR modeling
// //
//
 real gompertzMean(real time, real Asym, real offset, real growth){
    return(Asym*exp(-offset*growth^time));
 }
// will need to wait until stan rolls out vectorized support
// not sure if this really helps indexing the input vector for groups/time
 // vector gompertzMeanVec(real time, real Asym, real offset, real growth){
 //   int ntime = num_elements(time);
 //    return(rep_vector(Asym, ntime) .*
 //      exp(-rep_vector(offset, ntime) .*
 //          exp(time .* log(rep_vector(growth, ntime)))
 //      )
 //    );
 // }

 real fourParamLogistic(real time, real min, real max, real inflect_point, real slope){
    return(max + (min - max / (1 + (time / inflect_point ) ^ slope) ) );
 }

 vector growthMean(int n_data, int[] groups, real[] times, matrix theta, int parameterization_type){
  vector[n_data] out;

  if(parameterization_type == 0){
     out = UserDefinedGrowth(n_data, groups, times, theta);
     return(out);
   }

  for(idx in 1:n_data){
   if(parameterization_type == 1){
     if(cols(theta) != 3){
       reject("This growth parameterization requires 3 params");
     }
     out[idx] = gompertzMean(times[idx],
     theta[groups[idx], 1],
     theta[groups[idx], 2],
     theta[groups[idx], 3]);
   }else if(parameterization_type == 2){
     if(cols(theta) != 4){
       reject("The growth parameterization requires 4 params");
     }
     out[idx] = fourParamLogistic(times[idx],
     theta[groups[idx], 1],
     theta[groups[idx], 2],
     theta[groups[idx], 3],
     theta[groups[idx], 4]);
   }else{
     reject("Parameterization type is not valid", parameterization_type);
   }
  }

  return(out);
 }
