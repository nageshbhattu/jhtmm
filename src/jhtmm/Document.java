/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jhtmm;

import cc.mallet.types.Alphabet;
import cc.mallet.types.AlphabetCarrying;
import cc.mallet.types.FeatureSequence;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author user
 */
public class Document implements Serializable, AlphabetCarrying {
    Alphabet alphabet;
    ArrayList<FeatureSequence> sentences;
    Document(Alphabet dict){
        alphabet = dict;
        sentences = new ArrayList<>();
    }
    
    void addSequence(FeatureSequence fetSequence){
        sentences.add(fetSequence);
    }
    
    void addSequenceAt(FeatureSequence fetSequence,int index){
        sentences.add(index,fetSequence);
    }
    @Override
    public Alphabet getAlphabet() {
        return alphabet;
    }
    
    @Override
    public Alphabet[] getAlphabets() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
