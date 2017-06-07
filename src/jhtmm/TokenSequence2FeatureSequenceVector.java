/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jhtmm;

import cc.mallet.pipe.Pipe;
import cc.mallet.types.Alphabet;
import cc.mallet.types.FeatureSequence;
import cc.mallet.types.Instance;
import cc.mallet.types.TokenSequence;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author user
 */
public class TokenSequence2FeatureSequenceVector extends Pipe {
    public TokenSequence2FeatureSequenceVector (Alphabet dataDict)
	{
		super (dataDict, null);
	}

	public TokenSequence2FeatureSequenceVector ()
	{
		super(new Alphabet(), null);
	}
	
	public Instance pipe (Instance carrier)
	{
		TokenSequence ts = (TokenSequence) carrier.getData();
                Document fetSequenceList = 
                        new Document((Alphabet) getDataAlphabet());
                int preTi = 0;
                for (int ti = 0; ti < ts.size(); ti++) {
                    String token = ts.get(ti).getText();
                    Pattern p = Pattern.compile("([\\p{L}\\p{N}\\p{P}]+)[.;?!]|[.;?!]");
                    Matcher m = p.matcher(token);
                    if(!m.find()){
                        if(ti==ts.size()-1){
                            FeatureSequence ret = 
                            new FeatureSequence ((Alphabet)getDataAlphabet(), ti-preTi-1);
                            for (int fi = preTi; fi < ti; fi++) {
                                ret.add (ts.get(fi).getText());
                            }
                            fetSequenceList.addSequence(ret);
                        }
                        
                    }else {
                        FeatureSequence ret = 
                        new FeatureSequence ((Alphabet)getDataAlphabet(), ti-preTi-1);
                        for (int fi = preTi; fi < ti; fi++) {
                            ret.add (ts.get(fi).getText());
                        }
                        if(m.groupCount()!=0)
                            ret.add(m.group(1));
                        fetSequenceList.addSequence(ret);
                    } 
		}
                if(fetSequenceList.sentences.size()==0){
                    System.out.println("Document with zero ");
                }
		carrier.setData(fetSequenceList);
		return carrier;
	}
}
