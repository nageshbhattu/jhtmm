/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jhtmm;

import cc.mallet.types.Dirichlet;
import cc.mallet.types.FeatureSequence;
import cc.mallet.types.InstanceList;
import cc.mallet.util.Randoms;
import java.io.File;
/**
 *
 * @author user
 */
public class EM {
    // D -- numDocs
    // T -- numTopics
    // S -- numStates
    // V -- VocabularySize
    int numTopics;
    int numStates;
    int numDocs;
    int numTypes;
    double alpha;
    double beta;
    InstanceList docs;
    int numIterations;
    double logLikelihood;
    double[][] phi; // T X V 
    double[][] theta;       // estimated theta D X T
    double epsilon;                      // estimated epsilon
    double[][][] p_dwzpsi;  // The state probabilities, D X V X 
           // that is Pr(z,psi | d,w). The weird name reflects the arrangements
           // of these probabilities in this 3D array.
    void init(int topics, double alpha, double beta, int iters,  int seed,InstanceList ilist){
        numTopics = topics;
        this.alpha = alpha;
        this.beta = beta;
        numIterations = iters;
        docs = ilist;
        numDocs = ilist.size();
        numTypes = ilist.getDataAlphabet().size();
        RandInitParams(seed);
    };
    
    
  // EM inference for HTMM
    void Infer(){
        for (int iter = 0; iter < numIterations; iter++) {
          EStep();
          MStep();
          System.out.println("iteration " + iter + ", loglikelihood = " + logLikelihood );
        }
    }
    void EStep(){
        logLikelihood = 0;
        for (int di = 0; di < docs.size(); di++) {
          double ll = EStepSingleDoc(di);
          logLikelihood += ll;
        }
        // Here we take into account the priors of the parameters when computing
        // the likelihood:
        IncorporatePriorsIntoLikelihood();
    };
      // Computes expectations and contribution to the likelihood for a single
    // document.
    double EStepSingleDoc(int d){
        Document document = (Document) docs.get(d).getData();
        int numSentences = document.sentences.size();
        double[][] locals = new double[numSentences][numTopics];
        double local_ll = ComputeLocalProbsForDoc(document, locals);
        double[] init_probs = new double[numTopics*2];
        for (int ti = 0; ti < numTopics; ti++) {
            init_probs[ti] = theta[d][ti];
            init_probs[ti+numTopics] = 0;  // Document must begin with a topic transition.
        }
        FastRestrictedHMM f = new FastRestrictedHMM(epsilon,locals,theta[d],init_probs);
        
        f.ForwardBackward();
        f.ComputeLogLikelihood();
        return local_ll + f.getLogLikelihood();
    };

    // Compute the emission (local) probabilities for a certain document.
    double ComputeLocalProbsForDoc(Document doc, double[][] locals){
        double ll= 0.0;
        for (int sIndex = 0; sIndex < doc.sentences.size(); sIndex++) {
            

            // This method is used to compute local probabilities for a word or for a
            // sentence.
            // Actually, we compute potentials rather than emission probabilities
            // because we have to normalize them

            for (int ti = 0; ti < numTopics; ti++) {
              locals[sIndex][ti] = 1;
            }
            Normalize(numTopics, locals[sIndex]);
            ll += Math.log(numTopics);
            FeatureSequence sent = doc.sentences.get(sIndex);
            double min=0;
            double norm = 0;
            for (int wi = 0; wi < sent.size(); wi++) {
              int word = sent.getIndexAtPosition(wi);
              for (int ti = 0; ti < numTopics; ti++) {
                locals[sIndex][ti] += Math.log(phi[ti][word]);
                if(locals[sIndex][ti]<min){
                    min = locals[sIndex][ti];
                }
              }
            }
            for(int ti = 0;ti<numTopics;ti++){
                locals[sIndex][ti] -= min;
                locals[sIndex][ti] = Math.exp(locals[sIndex][ti]);
                norm += locals[sIndex][ti];
            }    
            Normalize(norm, locals[sIndex]);  // to prevent underflow
            ll += Math.log(norm);
        }
        return ll;
    }
    
    void Normalize(double factor, double[] vec){
        for(int vi = 0;vi<vec.length;vi++){
            vec[vi] = vec[vi]/factor;
        }
    };

  // Saves all parameters and the distribution on hidden states.
    void SaveAll(String base_name){};


  // Methods:
  // E-step of the algorithm
  


  // M-step of the algorithm
    void MStep(){
        FindEpsilon();
        FindPhi();
        FindTheta();
    };

  // Finds the MAP estimator for all thetas.
    void FindTheta(){
      for (int di = 0; di < docs.size(); di++) {
        FindSingleTheta(di);
      }
    };

  // Finds the MAP estimator for theta_d.
    void FindSingleTheta(int di){
        double norm = 0;
        double[] Cdz = new double[numTopics];
        CountTopicsInDoc(di, Cdz);
        for (int ti = 0; ti < numTopics; ti++) {
          theta[di][ti] = Cdz[ti] + alpha - 1;
          norm += theta[di][ti];
        }
        Normalize(norm, theta[di]);
    };

    
  // Finds the MAP estimator for epsilon.
    void FindEpsilon(){
        int total = 0;
        double lot = 0;
        for (int di = 0; di < docs.size(); di++) {
        //  we start counting from the second item in the document
            Document document = (Document) docs.get(di).getData();
            int numSentences= document.sentences.size();
            for (int sIndex = 1; sIndex < numSentences; sIndex++) {
              for (int ti = 0; ti < numTopics; ti++) {
                // only psi=1
                lot += p_dwzpsi[di][sIndex][ti];
              }
            }
            total += numSentences-1;      // Psi is always 1 for the first
                                              // word/sentence
        }
        epsilon = lot/total;  
    };

  // Finds the MAP estimator for all phi.
    void FindPhi(){
        double[][] Czw  = new double[numTopics][numTypes];
        CountTopicWord(Czw);   // Czw is allocated and initialized to 0
        for (int ti = 0; ti < numTopics; ti++) {
          double norm = 0;
          for (int wi = 0; wi < numTypes; wi++) {
            phi[ti][wi] = Czw[ti][wi] + beta - 1;
            norm += phi[ti][wi];
          }
          Normalize(norm, phi[ti]);
        }
    };

    // Counts (the expectation of) how many times a topic was drawn from
    // theta in a certain document.
    void CountTopicsInDoc(int di, double[] Cdz){
        Document document = (Document) docs.get(di).getData();
        int numSenteces = document.sentences.size();
        for (int sIndex = 0; sIndex < numSenteces; sIndex++) {
            for (int ti = 0; ti < numTopics; ti++) {
          // only psi=1
               Cdz[ti]+= p_dwzpsi[di][sIndex][ti];
            }
        }
    };

  // Counts (the expectation of) how many times the each pair (topic, word)
  // appears in a certain sentence.
    void CountTopicWordInSentence(FeatureSequence sent,
                                double[] topic_probs,
                                double[][] Czw){
    // Iterate over all the words in a sentence
        for (int wi = 0; wi < sent.size(); wi++) {
          int word = sent.getIndexAtPosition(wi);
          for (int ti = 0; ti < numTopics; ti++) {
            // both psi=1 and psi=0
            Czw[ti][wi] += topic_probs[ti]+topic_probs[ti+numTopics];
          }
        }  
    };

  // Counts (the expectation of) how many times the each pair (topic, word)
  // appears in the whole corpus.
    void CountTopicWord(double[][] Czw){
        for (int di = 0; di < docs.size(); di++) {
        Document document = (Document) docs.get(di).getData();
        int numSenteces = document.sentences.size();
        for (int sIndex = 0; sIndex < numSenteces; sIndex++) {
              CountTopicWordInSentence(document.sentences.get(sIndex), 
                      p_dwzpsi[di][sIndex], Czw);
            }
        }
    };
 
  // Saves theta in a file with extension .theta
    void SaveTheta(String fname){


    };

  // Saves phi in a file with extension .phi.
    void SavePhi(String fname){

    };

  // Saves epsilon in a file with extension .eps.
    void SaveEpsilon(String fname){

    };

  // Saves the latent states probabilities.
    void SaveTopicTransProbs(String fname){

    };

  // Saves the log likelihood.
    void SaveLogLikelihood(String fname){

    };

    // Randomly initializes the parameters (phi, theta and epsilon) and allocates
    // space for the state probabilities.
    // See comments above (for init) regarding the random seed.
    
    void RandInitParams(int seed) {
        // If seed is 0 use the clock to initialize the random generator
        
        Randoms random = new Randoms(seed);
        epsilon = random.nextDouble();
        theta = new double[numDocs][];

        Dirichlet dDirichlet = new Dirichlet(numTopics);
        for (int di = 0; di < docs.size(); di++) {
          theta[di] = dDirichlet.nextDistribution();
        }

        phi = new double[numTypes][];
        Dirichlet tDirichlet = new Dirichlet(numTypes);
        for(int ti =0;ti<numTopics;ti++){
            phi[ti] = tDirichlet.nextDistribution();
        }

        p_dwzpsi = new double[docs.size()][][];
        for (int di = 0; di < docs.size(); di++) {
            Document document = (Document) docs.get(di).getData();
            int numSentences= document.sentences.size();
            if(numSentences == 0){
                System.out.println("No Sentences in the document "+ di);
            }
            p_dwzpsi[di] = new double[numSentences][numTopics*2];
            
        }
    }


    
  // Given the parameters, finds the most probable sequence of hidden states.
    

  // Adds the Dirichlet priors to the likelihood computation.
    void IncorporatePriorsIntoLikelihood(){
        for (int di = 0; di < docs.size(); di++) {
            for (int ti = 0; ti < numTopics; ti++) {
              logLikelihood += (alpha-1)*Math.log(theta[di][ti]);
            }
          }
          // The prior on phi, assuming a symmetric Dirichlet distirubiton
          for (int ti = 0; ti < numTopics; ti++) {
            for (int wi = 0; wi < numTypes; wi++) {
              logLikelihood += (beta-1)*Math.log(phi[ti][wi]);
            }
          }
    };

  // I'm not sure if that's the place for this method!
  // Data members:
    
    public static void main(String[] args){
        InstanceList ilist = InstanceList.load(new File(args[0]));
        EM em = new EM();
        //void init(int topics, double alpha, double beta, int iters,  
        // int seed,InstanceList ilist)
        int numTopics = 10;
        double alpha = 1 + 50/numTopics;
        double beta = 1.01;
        int numIterations = 100;
        int seed = 1111;
        em.init(numTopics,alpha,beta,numIterations,seed,ilist);
        em.Infer();
    }

}
