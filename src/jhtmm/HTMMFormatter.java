/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jhtmm;

import cc.mallet.types.FeatureSequence;
import cc.mallet.types.InstanceList;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author user
 */
public class HTMMFormatter {
    public static void main(String[] args) throws IOException{
        InstanceList ilist = InstanceList.load(new File(args[0]));
        System.out.println("Size of the Vocabulary " + ilist.getDataAlphabet().size());
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File("filenames.txt")));
        for(int di = 0;di<ilist.size();di++){
            bw.write("input/"+di+"\n");
            BufferedWriter bwDoc = new BufferedWriter(new FileWriter(new File("input/"+di)));
            
            Document d = (Document) ilist.get(di).getData();
            
            for(int sIndex = 0;sIndex<d.sentences.size();sIndex++){
                FeatureSequence fs = d.sentences.get(sIndex);
                bwDoc.write(fs.size() + " ");
                int[] features = fs.getFeatures();
                for(int wi = 0;wi<fs.size();wi++){
                    bwDoc.write(features[wi]+" ");
                    System.out.println("Htmm prob");
                }
                bwDoc.write("\n");
            }
            bwDoc.close();
        }
        bw.close();
    }
}
