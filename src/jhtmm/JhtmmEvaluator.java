/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jhtmm;

import cc.mallet.types.Dirichlet;
import cc.mallet.types.FeatureSequence;
import cc.mallet.types.Instance;
import cc.mallet.types.InstanceList;
import cc.mallet.util.Randoms;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author nageshbhattu
 */
public class JhtmmEvaluator {
    protected int numTopics; // Number of topics to be fit

	protected double[] alpha;	 // Dirichlet(alpha,alpha,...) is the distribution over topics
	protected double alphaSum;
	protected double beta;   // Prior on per-topic multinomial distribution over words
	protected double betaSum;
	

	protected double[][] phi; // indexed by <feature index, topic index>
	protected double[] tokensPerTopic; // indexed by <topic index>

	protected boolean printWordProbabilities = false;
        protected int numIterations =100;
        protected double logLikelihood = 0.0;
	protected int wordStopIndex = 4;
        protected double[][] theta;
        protected double epsilon;
        protected int lastSentIndex;
        protected int lastWordIndex;
        protected double[][] sprobs;    
        int numStates;
    private double[][] fprobs;
    private double[][] bprobs;
    private double[] norm_factor;
    int testDocSize;
	public JhtmmEvaluator (int numTopics,
                              double[] alpha, double alphaSum,
                              double beta,
                              double[][] phi, 
                              int wordStopIndex,
                              double epsilon) {

            this.numTopics = numTopics;
            this.numStates = numTopics*2;

            this.phi = phi;
            this.alphaSum = alphaSum;
            this.alpha = alpha;
            this.beta = beta;
            this.betaSum = beta * phi.length;
            this.wordStopIndex = wordStopIndex;
            this.epsilon = epsilon;
	}
        
        int testSentIndex = 0;

    JhtmmEvaluator() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

     
        public void evaluateTM(InstanceList ilist){
            Dirichlet dDirichlet = new Dirichlet(numTopics);
            theta = new double[ilist.size()][numTopics];
            int dIndex = 0;
            for(Instance inst:ilist){
                theta[dIndex] = dDirichlet.nextDistribution();
                //calculate theta_new using EM algorithm
                for (int iter = 0; iter < numIterations; iter++) {
                    EStepSingleDoc(inst,dIndex);
                    
                    
                //
                
                }
                dIndex++;
            }
        }
      // Computes expectations and contribution to the likelihood for a single
    // document.
        
        double [][] locals;
        void splitTestSentence(Document document){
            int stopIndex = wordStopIndex;
            for(int sIndex = 0;sIndex< document.sentences.size() && stopIndex>0;sIndex++){
                int sentenceLength = document.sentences.get(sIndex).size();
                if(stopIndex>=document.sentences.get(sIndex).size())
                    stopIndex-= sentenceLength;
                else{
                    
                    FeatureSequence fetSequence = document.sentences.get(sIndex);
                    int[]features =  fetSequence.getFeatures();
                    
                    FeatureSequence leftFeatureSequence = new FeatureSequence(fetSequence.getAlphabet(), Arrays.copyOfRange(features, 0, stopIndex));
                    FeatureSequence rightFeatureSequence = new FeatureSequence(fetSequence.getAlphabet(), Arrays.copyOfRange(features, stopIndex, sentenceLength));
                    document.sentences.set(sIndex, leftFeatureSequence);
                    document.sentences.add(sIndex+1,rightFeatureSequence);
                    testSentIndex = sIndex;
                    sIndex++;
                    break;
                }
            }
            int docSize =0;
            for(int sIndex = 0;sIndex<document.sentences.size();sIndex++){
                docSize += document.sentences.get(sIndex).getLength();
            }
            testDocSize = docSize-wordStopIndex;
        }
        
        void EStepSingleDoc(Instance instance,int dIndex){
            Document document = (Document) instance.getData();
            splitTestSentence(document);
            
            int numSentences = document.sentences.size();
            
            locals = new double[numSentences][numTopics];
            
            ComputeLocalProbsForDoc(document, locals);
            
            double[] init_probs = new double[numTopics*2];
            
            for (int ti = 0; ti < numTopics; ti++) {
                init_probs[ti] = theta[dIndex][ti];
                init_probs[ti+numTopics] = 0;  // Document must begin with a topic transition.
            }
            
            FastRestrictedHMMTest f = new FastRestrictedHMMTest(epsilon,locals,theta[dIndex],init_probs,lastSentIndex,testDocSize);

             // Perform M Step of EM algorithm for the computation of theta_d
             double[][] sprobs = f.sprobs;
             for(int ti = 0;ti<numTopics;ti++){
                 theta[dIndex][ti] = 0.0;
             }

             for(int sIndex = 0;sIndex<=lastSentIndex;sIndex++){
                 for(int ti = 0;ti< numTopics;ti++){
                     theta[dIndex][ti]+=sprobs[sIndex][ti];
                 }
             }
             // Use the theta_d computed in the previous step for the remaining part of the document to compute the perplexity
             

        };
    void ComputeLocalProbsForDoc(Document doc, double[][] locals){
        for (int sIndex = 0; sIndex < doc.sentences.size(); sIndex++) {
            

            // This method is used to compute local probabilities for a word or for a
            // sentence.
            // Actually, we compute potentials rather than emission probabilities
            // because we have to normalize them

            for (int ti = 0; ti < numTopics; ti++) {
              locals[sIndex][ti] = 1;
            }
            Normalize(numTopics, locals[sIndex]);
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
        }
    }
    
    
    void Normalize(double factor, double[] vec){
        for(int vi = 0;vi<vec.length;vi++){
            vec[vi] = vec[vi]/factor;
        }
    };
     public static void main(String[] args){
         
        InstanceList ilist = InstanceList.load(new File(args[0]));
        double [] proportions = {0.75,0.25};
        InstanceList[] ilists = ilist.split(proportions);
        EM em = new EM();
        //void init(int topics, double alpha, double beta, int iters,  
        // int seed,InstanceList ilist)
        int numTopics = 10;
        double alpha = 1 + 50/numTopics;
        double beta = 1.01;
        int numIterations = 100;
        int seed = 1111;
        em.init(numTopics,alpha,beta,numIterations,seed,ilists[0]);
        em.Infer();
        double[] alphaArr = new double[numTopics];
        for(int ti = 0;ti<numTopics;ti++)
            alphaArr[ti] = alpha;
        
        int wordStopIndex = 64;
        double epsilon = 0.3762;
        JhtmmEvaluator jhtmmevaluator = new JhtmmEvaluator(numTopics, alphaArr, alpha*numTopics,beta,em.phi, wordStopIndex, epsilon);
        //void init(int topics, double alpha, double beta, int iters,  
        // int seed,InstanceList ilist)
        jhtmmevaluator.evaluateTM(ilists[1]);
    }

    
    
    
    };
    

