/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jhtmm;

import java.util.ArrayList;

/**
 *
 * @author user
 */
public class FastRestrictedHMMTest {
    int numTopics;
    int numStates;
    int numSentences;
    double epsilon;
    double[][] local; // numSentences * numTopics
    double[] theta; // numTopics
    double[] pi; // numStates
    double[][] sprobs; // numSentences * numStates
    double[][] alpha; // numSentences * numStates forward probabilities
    double[][] beta ; // numSentences * numStates backward probabilities
    double[]norm_factor; // numSentences 
    int begSentIndex;
    int endSentIndex;
    int begWordIndex;
    int endWordIndex;
    int testDocSize;
    
    
    public FastRestrictedHMMTest(double epsilon,double[][]local,double[] theta,double[]pi,
                                 int sentIndex, int testDocSize){
        this.pi = pi;
        this.theta = theta;
        this.epsilon = epsilon;
        this.local = local;
        numTopics = theta.length;
        numStates = numTopics*2;
        numSentences = local.length;
        if(numSentences==0){
            System.out.println("Something Wrong");
        }
        alpha = new double[numSentences][numStates];
        beta = new double[numSentences][numStates];
        norm_factor = new double[numSentences];
        sprobs = new double[numSentences][numStates];
        begSentIndex = 0;
        endSentIndex = sentIndex;
        this.testDocSize= testDocSize; 
    }
    public double computePerplexity(){
        InitAlpha();
        
        ComputeAlphas(1,endSentIndex);
        
        InitBeta();
        
        ComputeBetas();
        
        CombineAllProbs();
        
        ComputeAlphas(endSentIndex+1,numSentences-1);
        
        double perplexity = 0.0;
        for(int ti = 0;ti<numStates;ti++){
            perplexity += alpha[endSentIndex][ti];
        }
        
        perplexity = Math.exp(-Math.log(perplexity)/testDocSize);
        
        return perplexity;
    }
    
    
    public double[][] getStateProbabilities(){
        return sprobs;
    }
    
    void CombineAllProbs(){
        for (int sIndex = begSentIndex; sIndex <= endSentIndex; sIndex++) { // sIndex is sentence Index
            double norm = 0;
            for (int si = 0; si < numStates; si++) { //si is state Index
                sprobs[sIndex][si] = alpha[sIndex][si]*beta[sIndex][si];
                norm += sprobs[sIndex][si];
            }
            Normalize(norm,sprobs[sIndex]);
        }
    };
    
    
  // This method normalizes the vector vec to according to the factor norm.
    void Normalize(double factor, double[] vec){
        for(int vi = 0;vi<vec.length;vi++){
            vec[vi] = vec[vi]/factor;
        }
    };

    // This method initializes alpha[0] to be Pr(y_0 | q_0; param) * Pr(q_0)
    // norm is the normalizing factor used.
    // The other parameters are as defined above.
    void InitAlpha()  {
        double normalizer = 0.0;
        normalizer = 0;
        for (int i = 0; i < numTopics; i++) {
            alpha[begSentIndex][i] = local[0][i]*pi[i];
            alpha[begSentIndex][i+numTopics] = local[0][i]*pi[i+numTopics];
            normalizer += alpha[begSentIndex][i] + alpha[begSentIndex][i+numTopics];
        }
        Normalize(normalizer, alpha[begSentIndex]);
        norm_factor[0] = normalizer;
    };

  // This method initializes beta[T-1] (the argument beta_T_1) to be all ones
  // and then normalizes it according to the normalization of alpha[T-1].
  // norm is the normalizing factor used.
  // The other parameters are as defined above.
  void InitBeta() {
      for(int ti = 0;ti<numTopics;ti++){
          beta[endSentIndex][ti] = 1;
      }
      Normalize(norm_factor[endSentIndex],beta[endSentIndex]);
  };

  // This method computes alpha[i] for i>=1 (after alpha[0] has been
  // initialized).
  // norm_factor is the ArrayList of normalizing factors used along the chain.
  // The other parameters are defined above.
  void ComputeAlphas(int begSIndex, int endSIndex){
        double normalizer = 0.0;
        
        for (int sIndex = begSIndex; sIndex < endSIndex; sIndex++) {
            normalizer = 0.0;
            
            for (int si = 0; si < numTopics; si++) {
                alpha[sIndex][si] = epsilon*theta[si]*local[sIndex][si];  // regardless of the previous
                   // topic - remember that sum_k alpha[t-1][k] is 1 (because of the norm).
                alpha[sIndex][si+numTopics] = (1-epsilon)*
                        (alpha[sIndex-1][si] + alpha[sIndex-1][si+numTopics])*
                        local[sIndex][si];
                norm_factor[sIndex] += alpha[sIndex][si]+alpha[sIndex][si+numTopics];
            }
            if(norm_factor[sIndex]==0){
                System.out.println("wrong sIndex");
            }
            Normalize(norm_factor[sIndex], alpha[sIndex]);
        }
  };


  // Computes beta[i] for i<=T-2 (beta[T-1] has been initialized). Use same
  // normalizing factors as in the alpha computation.
  // norm_factor is the vector of normalizing factors used along the chain.
  // The other parameters are as defined above.
  void ComputeBetas(){
    for (int sIndex = endSentIndex-1; sIndex >= begSentIndex; sIndex--) {
        double trans_sum = 0;
        for (int ti = 0; ti < numTopics; ti++) {
            trans_sum += epsilon*theta[ti]*local[sIndex+1][ti]*beta[sIndex+1][ti];
        }

        for (int ti = 0; ti < numTopics; ti++) {
        // Recall that beta_t1[s] == beta_t1[s+topics_]
            beta[sIndex][ti] = trans_sum + (1-epsilon)*local[sIndex+1][ti]*beta[sIndex+1][ti];
            beta[sIndex][ti]/=norm_factor[sIndex];
            beta[sIndex][ti+numTopics] = beta[sIndex][ti];
        }
    }
  };
}
