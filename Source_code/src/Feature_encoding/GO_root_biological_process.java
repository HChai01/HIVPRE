package Feature_encoding;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.Writer;

/*
 * This Java code is designed to map the Gene Ontology term to the child term of biological process (GO:0008150) through the derivation tree.
 * 
 * Feature list:
 * F1: GO:0000003	reproduction
 * F2: GO:0001906	cell_killing
 * F3: GO:0002376	immune_system_process
 * F4: GO:0006791	sulfur_utilization
 * F5: GO:0006794	phosphorus_utilization
 * F6: GO:0007610	behavior
 * F7: GO:0008152	metabolic_process
 * F8: GO:0008283	cell_population_proliferation
 * F9: GO:0009758	carbohydrate_utilization
 * F10: GO:0009987	cellular_process
 * F11: GO:0015976	carbon_utilization
 * F12: GO:0019740	nitrogen_utilization
 * F13: GO:0022414	reproductive_process
 * F14: GO:0022610	biological_adhesion
 * F15: GO:0023052	signaling
 * F16: GO:0032501	multicellular_organismal_process
 * F17: GO:0032502	developmental_process
 * F18: GO:0040007	growth
 * F19: GO:0040011	locomotion
 * F20: GO:0043473	pigmentation
 * F21: GO:0044848	biological_phase
 * F22: GO:0048511	rhythmic_process
 * F23: GO:0050896	response_to_stimulus
 * F24: GO:0051179	localization
 * F25: GO:0051704	multi-organism_process
 * F26: GO:0065007	biological_regulation
 * F27: GO:0071840	cellular_component_organization_or_biogenesis
 * F28: GO:0098743	cell_aggregation
 * F29: GO:0098754	detoxification
 * 
 * Citation: Chai H, Gu Q, Hughes J, Robertson DL (2022) In silico prediction of HIV-1-host molecular interactions and their directionality. PLoS Comput Biol 18(2): e1009720. https://doi.org/10.1371/journal.pcbi.1009720.
 * 
 * */

public class GO_root_biological_process {
	// Path of the main folder
	static String topPath = "/Users/imac/Documents/HIVPRE/data/";
	
	static String sourcePath = topPath + "Gene_Ontology_example.txt";
	static String mapPath = topPath + "GO_BProot_pure.txt";
	static String outPath = topPath + "GO_root_biological_process.txt";
	static int rootnum = 29;
	
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
