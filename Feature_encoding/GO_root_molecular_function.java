package Feature_encoding;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.Writer;

/*
 * This Java code is designed to map the Gene Ontology term to the child term of molecular function (GO:0003674) through the derivation tree.
 * 
 * Feature list:
 * F1:	GO:0003824	catalytic_activity
 * F2:	GO:0005198	structural_molecule_activity
 * F3:	GO:0005215	transporter_activity
 * F4:	GO:0005488	binding
 * F5:	GO:0016209	antioxidant_activity
 * F6:	GO:0031386	protein_tag
 * F7:	GO:0038024	cargo_receptor_activity
 * F8:	GO:0044183	protein_folding_chaperone
 * F9:	GO:0045182	translation_regulator_activity
 * F10:	GO:0045735	nutrient_reservoir_activity
 * F11:	GO:0060089	molecular_transducer_activity
 * F12:	GO:0090729	toxin_activity
 * F13:	GO:0098772	molecular_function_regulator
 * F14:	GO:0140104	molecular_carrier_activity
 * F15:	GO:0140110	transcription_regulator_activity
 * F16:	GO:0140299	small_molecule_sensor_activity
 * 
 * Citation: Haiting Chai, Quan Gu, Joseph Hughes and David L. Robertson (202X), In silico prediction of HIV-1-host molecular interactions and their directionality, PLOS Computational Biology, XX(X): XXX-XXX, PMID: XXXXXXXX.
 * 
 * */

public class GO_root_molecular_function {
	// Path of the main folder
	static String topPath = "/Users/imac/Documents/HIVPRE/data/";
	
	static String sourcePath = topPath + "Gene_Ontology_example.txt";
	static String mapPath = topPath + "GO_MFroot_pure.txt";
	static String outPath = topPath + "GO_root_molecular_function.txt";
	static int rootnum = 16;
	
	static int readpointer;
	static int temporarypointer;
	static int rootpointer;
	static int respointer;
	
	static String Info[][] = new String[100000][3];
	static int res[] = new int[100];
	
	public static void readmap(String Path){
		String[] Recordcut;
		readpointer = 1;
		try {
			String encoding = "UTF-8";
            File listfile=new File(Path);
            if(listfile.isFile() && listfile.exists()){ 
            	InputStreamReader read = new InputStreamReader(
                new FileInputStream(listfile),encoding);                    
                BufferedReader bufferedReader = new BufferedReader(read);
                String lineTxt = "N/A";
                lineTxt = bufferedReader.readLine();
                lineTxt = bufferedReader.readLine();
                lineTxt = bufferedReader.readLine();
                lineTxt = bufferedReader.readLine();
                while((lineTxt = bufferedReader.readLine()) != null){
                		Recordcut = lineTxt.split("\\s{1,}");
                		Info[readpointer][0] = Recordcut[0];
                		Info[readpointer][1] = Recordcut[1];
                		Info[readpointer][2] = Recordcut[2];
                		readpointer = readpointer + 1;
                }
                Info[0][0] = readpointer - 1 + "";
                write(outPath,"Index\tHeadinfo");
                for(readpointer=1;readpointer<rootnum+1;readpointer++) {
                		write(outPath,"\t"+Info[readpointer][1]);
                }
                write(outPath,"\n");
                read.close();
            }else{
            	System.out.println("Source File missing!");
            }
        } catch (Exception e) {
            System.out.println("Source Read error!");
            e.printStackTrace();
        }
    }
	
	public static void readsource(String Path){
		String[] Recordcut;
		readpointer = 0;
		try {
			String encoding = "UTF-8";
            File listfile=new File(Path);
            if(listfile.isFile() && listfile.exists()){ 
            	InputStreamReader read = new InputStreamReader(
                new FileInputStream(listfile),encoding);                    
                BufferedReader bufferedReader = new BufferedReader(read);
                String lineTxt = "N/A";
                while((lineTxt = bufferedReader.readLine()) != null){
            			readpointer = readpointer + 1;
                		Recordcut = lineTxt.split("\\s{1,}");
                		write(outPath,"G"+readpointer);
                		if(Recordcut.length!=2) {
                			System.out.println("Format unsupported, please check examples in 'Gene_Ontology_example.txt'.");
                			write(outPath,"\tFormat unsupported.\n");
                			continue;
                		}
                		write(outPath,"\t"+Recordcut[0]);
                		mapGOroot(Recordcut[0],Recordcut[1]);              		              		
                }
                read.close();
            }else{
            	System.out.println("Source File missing!");
            }
        } catch (Exception e) {
            System.out.println("Source Read error!");
            e.printStackTrace();
        }
    }
	
	public static void mapGOroot(String headinfo,String clue){
		String[] cluecut;
		for(respointer=1;respointer<rootnum+1;respointer++) {
			res[respointer] = 0;
	    }
		if(clue.contains(";")&&clue.contains(",")) {
			System.out.println("Septation chaos, please check content in: " + headinfo);
			write(outPath,"\tSeptation chaos.\n");
		}else {
			if(clue.contains(";")) {
				cluecut = clue.split(";");
				for(temporarypointer=0;temporarypointer<cluecut.length;temporarypointer++) {
					checkroot(cluecut[temporarypointer]);					
				}
				for(respointer=1;respointer<rootnum+1;respointer++) {
					write(outPath,"\t" + res[respointer] + "");
				}
				write(outPath,"\n");
			}else {
				if(clue.contains(",")) {
					cluecut = clue.split(",");
					for(temporarypointer=0;temporarypointer<cluecut.length;temporarypointer++) {
						checkroot(cluecut[temporarypointer]);						
					}
					for(respointer=1;respointer<rootnum+1;respointer++) {
						write(outPath,"\t" + res[respointer] + "");
					}
					write(outPath,"\n");
				}else {
					System.out.println("Unrecognized septation, please check content in: " + headinfo);
					write(outPath,"\tSeptation unrecognized.\n");
				}
			}
		}	
	}	
	
	public static void checkroot(String GOterm){
		for(rootpointer=1;rootpointer<Integer.parseInt(Info[0][0])+1;rootpointer++) {
			if(GOterm.equals(Info[rootpointer][0])) {
				for(respointer=1;respointer<rootnum+1;respointer++) {
		        		if(Info[rootpointer][2].contains(Info[respointer][0])) {
		        			res[respointer] = 1;
		        		}
		        }
				break;
			}
		}
	}
	
	public static void write(String Path, String textstr){
		try {
			File textfile = new File(Path);	
			Writer output = null;
			output = new FileWriter(textfile,true);
			output.write(textstr);
			output.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public static void main(String argv[]){ 
		readmap(mapPath);
		readsource(sourcePath);
	}
}
