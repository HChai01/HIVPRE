package Feature_encoding;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.DecimalFormat;

/*
 * This Java code is designed to calculate the codon usage in the coding sequence.
 * 
 * Feature list:
 * F1: TTT_(F)
 * F2: TTC_(F)
 * F3: TTA_(L)
 * F4: TTG_(L)
 * F5: CTT_(L)
 * F6: CTC_(L)
 * F7: CTA_(L)
 * F8: CTG_(L)
 * F9: ATT_(I)
 * F10: ATC_(I)
 * F11: ATA_(I)
 * F12: ATG_(M)
 * F13: GTT_(V)
 * F14: GTC_(V)
 * F15: GTA_(V)
 * F16: GTG_(V)
 * F17: TCT_(S)
 * F18: TCC_(S)
 * F19: TCA_(S)
 * F20: TCG_(S)
 * F21: AGT_(S)
 * F22: AGC_(S)
 * F23: CCT_(P)
 * F24: CCC_(P)
 * F25: CCA_(P)
 * F26: CCG_(P)
 * F27: ACT_(T)
 * F28: ACC_(T)
 * F29: ACA_(T)
 * F30: ACG_(T)
 * F31: TAT_(Y)
 * F32: TAC_(Y)
 * F33: GCT_(A)
 * F34: GCC_(A)
 * F35: GCA_(A)
 * F36: GCG_(A)
 * F37: CAT_(H)
 * F38: CAC_(H)
 * F39: CAA_(Q)
 * F40: CAG_(Q)
 * F41: AAT_(N)
 * F42: AAC_(N)
 * F43: AAA_(K)
 * F44: AAG_(K)
 * F45: GAT_(D)
 * F46: GAC_(D)
 * F47: GAA_(E)
 * F48: GAG_(E)
 * F49: TGT_(C)
 * F50: TGC_(C)
 * F51: TGG_(W)
 * F52: CGT_(R)
 * F53: CGC_(R)
 * F54: CGA_(R)
 * F55: CGG_(R)
 * F56: AGA_(R)
 * F57: AGG_(R)
 * F58: GGT_(G)
 * F59: GGC_(G)
 * F60: GGA_(G)
 * F61: GGG_(G)
 * F62: TAA_(N/A)
 * F63: TAG_(N/A)
 * F64: TGA_(N/A)
 * 
 * Citation: Haiting Chai, Quan Gu, Joseph Hughes and David L. Robertson (202X), In silico prediction of HIV-1-host molecular interactions and their directionality, PLOS Computational Biology, XX(X): XXX-XXX, PMID: XXXXXXXX.
 * 
 * */

public class Codon_usage_CDS {
	// Path of the main folder
	static String topPath = "/Users/imac/Documents/HIVPRE/data/";
	
	static String CDSPath = topPath + "DNA_CDs_example.txt";
	static String outPath = topPath + "Condon_usage(ver.CDS).txt";
	
	static int readpointer;
	static int temporarypointer;
	static String strtxt = "";
	static String seq = "";
	static int[] codonINFO = new int[65];
	
	static DecimalFormat decimalFormat = new DecimalFormat("0.0000");
	
	public static void readCDS(String Path){
		String headinfo;
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
                if(!lineTxt.contains(">")) {
                		System.out.println("Please check the headline of your CDS.");
                }else {
                		headinfo = lineTxt.trim().replace(">", "");;
                		strtxt = "S" + readpointer + "\t" + headinfo;
                		seq = "";
                		while((lineTxt = bufferedReader.readLine()) != null){
                			if(lineTxt.contains(">")) {
                				write(outPath,strtxt);
                				calculatecodon(seq,headinfo);              				
                				readpointer = readpointer + 1;
                				headinfo = lineTxt.trim().replace(">", "");
                				strtxt = "S" + readpointer + "\t" + headinfo;
                        		seq = "";
                			}else {
                				seq = seq + lineTxt.trim();
                			}
                		}
                		write(outPath,strtxt);
                		calculatecodon(seq,headinfo);
                }
                read.close();
            }else{
            	System.out.println("CDS File missing!");
            }
        } catch (Exception e) {
            System.out.println("CDS File Read error!");
            e.printStackTrace();
        }
    }
	
	public static void calculatecodon(String query,String id){
		String stopcodon;
		String startcodon;
		String clue = "";
		int key = 0;
		for(temporarypointer=0;temporarypointer<65;temporarypointer++) {
			codonINFO[temporarypointer] = 0;
  		}
		char a[] = query.toCharArray();
		if(a.length%3!=0 || a.length <6) {
			System.out.println("CDS is incompleted, please check the content in: " + id);
			write(outPath,"\tCDS incompleted\n");
		}else {
			stopcodon = a[a.length-3] + "" + a[a.length-2] + "" + a[a.length-1] + "";
			if(!(stopcodon.equals("TAA")||stopcodon.equals("TAG")||stopcodon.equals("TGA"))) {
				System.out.println("Stop codon is missing, please check the content in: " + id);
				write(outPath,"\tStop codon missing\n");
			}else {
				startcodon = a[0] + "" + a[1] + "" + a[2] + "";
				if(!startcodon.equals("ATG")) {
					System.out.println("Start codon is abnormal, please check the content in: " + id);
					write(outPath,"\tStart codon abnormal");
				}
				for(temporarypointer=0;temporarypointer<a.length/3;temporarypointer++) {
	        			clue = a[temporarypointer*3] + "" + a[temporarypointer*3+1] + "" + a[temporarypointer*3+2] + "";
	        			countcodon(clue);
	        			if(temporarypointer < (a.length/3)-1 && (clue.equals("TAA")||clue.equals("TAG")||clue.equals("TGA"))) {
	        				System.out.println("Abnormal stop codons '" + clue + "' is detected, please check the content in: " + id);
	        				key = 1;
	        			}
        			}
				if(key == 1) {
					write(outPath,"\tMultiple stop codons detected");
				}else {
					write(outPath,"\tNormal");
				}
				
				for(temporarypointer=1;temporarypointer<65;temporarypointer++) {
	        			if(codonINFO[temporarypointer]==0) {
	        				write(outPath,"\t0.0000");
	        			}else {          			
	                		if(temporarypointer>=1 && temporarypointer<=2) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[1]+codonINFO[2])));
	                		}
	                		if(temporarypointer>=3 && temporarypointer<=8) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[3]+codonINFO[4]+codonINFO[5]+codonINFO[6]+codonINFO[7]+codonINFO[8])));
	                		}
	                		if(temporarypointer>=9 && temporarypointer<=11) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[9]+codonINFO[10]+codonINFO[11])));
	                		}
	                		if(temporarypointer==12) {
	                			write(outPath,"\t1.0000");
	                		}
	                		if(temporarypointer>=13 && temporarypointer<=16) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[13]+codonINFO[14]+codonINFO[15]+codonINFO[16])));
	                		}
	                		if(temporarypointer>=17 && temporarypointer<=22) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[17]+codonINFO[18]+codonINFO[19]+codonINFO[20]+codonINFO[21]+codonINFO[22])));
	                		}
	                		if(temporarypointer>=23 && temporarypointer<=26) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[23]+codonINFO[24]+codonINFO[25]+codonINFO[26])));
	                		}
	                		if(temporarypointer>=27 && temporarypointer<=30) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[27]+codonINFO[28]+codonINFO[29]+codonINFO[30])));
	                		}
	                		if(temporarypointer>=31 && temporarypointer<=32) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[31]+codonINFO[32])));
	                		}
	                		if(temporarypointer>=33 && temporarypointer<=36) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[33]+codonINFO[34]+codonINFO[35]+codonINFO[36])));
	                		}
	                		if(temporarypointer>=37 && temporarypointer<=38) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[37]+codonINFO[38])));
	                		}
	                		if(temporarypointer>=39 && temporarypointer<=40) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[39]+codonINFO[40])));
	                		}
	                		if(temporarypointer>=41 && temporarypointer<=42) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[41]+codonINFO[42])));
	                		}
	                		if(temporarypointer>=43 && temporarypointer<=44) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[43]+codonINFO[44])));
	                		}
	                		if(temporarypointer>=45 && temporarypointer<=46) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[45]+codonINFO[46])));
	                		}
	                		if(temporarypointer>=47 && temporarypointer<=48) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[47]+codonINFO[48])));
	                		}
	                		if(temporarypointer>=49 && temporarypointer<=50) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[49]+codonINFO[50])));
	                		}
	                		if(temporarypointer==51) {
	                			write(outPath,"\t1.0000");
	                		}
	                		if(temporarypointer>=52 && temporarypointer<=57) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[52]+codonINFO[53]+codonINFO[54]+codonINFO[55]+codonINFO[56]+codonINFO[57])));
	                		}
	                		if(temporarypointer>=58 && temporarypointer<=61) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[58]+codonINFO[59]+codonINFO[60]+codonINFO[61])));
	                		}
	                		if(temporarypointer>=62 && temporarypointer<=64) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[62]+codonINFO[63]+codonINFO[64])));
	                		}
	        			}
				}
				write(outPath,"\n");
			}
		}		
	}
	
	public static void countcodon(String codon) {
		switch(codon) {
			case	 "TTT":
				codonINFO[1] = codonINFO[1] + 1;
				break;
			case	 "TTC":
				codonINFO[2] = codonINFO[2] + 1;
				break;
			case	 "TTA":
				codonINFO[3] = codonINFO[3] + 1;
				break;
			case	 "TTG":
				codonINFO[4] = codonINFO[4] + 1;
				break;
			case	 "CTT":
				codonINFO[5] = codonINFO[5] + 1;
				break;
			case	 "CTC":
				codonINFO[6] = codonINFO[6] + 1;
				break;
			case	 "CTA":
				codonINFO[7] = codonINFO[7] + 1;
				break;
			case	 "CTG":
				codonINFO[8] = codonINFO[8] + 1;
				break;
			case	 "ATT":
				codonINFO[9] = codonINFO[9] + 1;
				break;
			case	 "ATC":
				codonINFO[10] = codonINFO[10] + 1;
				break;
			case	 "ATA":
				codonINFO[11] = codonINFO[11] + 1;
				break;
			case	 "ATG":
				codonINFO[12] = codonINFO[12] + 1;
				break;
			case	 "GTT":
				codonINFO[13] = codonINFO[13] + 1;
				break;
			case	 "GTC":
				codonINFO[14] = codonINFO[14] + 1;
				break;
			case	 "GTA":
				codonINFO[15] = codonINFO[15] + 1;
				break;
			case	 "GTG":
				codonINFO[16] = codonINFO[16] + 1;
				break;
			case	 "TCT":
				codonINFO[17] = codonINFO[17] + 1;
				break;
			case	 "TCC":
				codonINFO[18] = codonINFO[18] + 1;
				break;
			case	 "TCA":
				codonINFO[19] = codonINFO[19] + 1;
				break;
			case	 "TCG":
				codonINFO[20] = codonINFO[20] + 1;
				break;
			case	 "AGT":
				codonINFO[21] = codonINFO[21] + 1;
				break;
			case	 "AGC":
				codonINFO[22] = codonINFO[22] + 1;
				break;
			case	 "CCT":
				codonINFO[23] = codonINFO[23] + 1;
				break;
			case	 "CCC":
				codonINFO[24] = codonINFO[24] + 1;
				break;
			case	 "CCA":
				codonINFO[25] = codonINFO[25] + 1;
				break;
			case	 "CCG":
				codonINFO[26] = codonINFO[26] + 1;
				break;
			case	 "ACT":
				codonINFO[27] = codonINFO[27] + 1;
				break;
			case	 "ACC":
				codonINFO[28] = codonINFO[28] + 1;
				break;
			case	 "ACA":
				codonINFO[29] = codonINFO[29] + 1;
				break;
			case	 "ACG":
				codonINFO[30] = codonINFO[30] + 1;
				break;
			case	 "TAT":
				codonINFO[31] = codonINFO[31] + 1;
				break;
			case	 "TAC":
				codonINFO[32] = codonINFO[32] + 1;
				break;
			case	 "GCT":
				codonINFO[33] = codonINFO[33] + 1;
				break;
			case	 "GCC":
				codonINFO[34] = codonINFO[34] + 1;
				break;
			case	 "GCA":
				codonINFO[35] = codonINFO[35] + 1;
				break;
			case	 "GCG":
				codonINFO[36] = codonINFO[36] + 1;
				break;
			case	 "CAT":
				codonINFO[37] = codonINFO[37] + 1;
				break;
			case	 "CAC":
				codonINFO[38] = codonINFO[38] + 1;
				break;
			case	 "CAA":
				codonINFO[39] = codonINFO[39] + 1;
				break;
			case	 "CAG":
				codonINFO[40] = codonINFO[40] + 1;
				break;
			case	 "AAT":
				codonINFO[41] = codonINFO[41] + 1;
				break;
			case	 "AAC":
				codonINFO[42] = codonINFO[42] + 1;
				break;
			case	 "AAA":
				codonINFO[43] = codonINFO[43] + 1;
				break;
			case	 "AAG":
				codonINFO[44] = codonINFO[44] + 1;
				break;
			case	 "GAT":
				codonINFO[45] = codonINFO[45] + 1;
				break;
			case	 "GAC":
				codonINFO[46] = codonINFO[46] + 1;
				break;
			case	 "GAA":
				codonINFO[47] = codonINFO[47] + 1;
				break;
			case	 "GAG":
				codonINFO[48] = codonINFO[48] + 1;
				break;
			case	 "TGT":
				codonINFO[49] = codonINFO[49] + 1;
				break;
			case	 "TGC":
				codonINFO[50] = codonINFO[50] + 1;
				break;
			case	 "TGG":
				codonINFO[51] = codonINFO[51] + 1;
				break;
			case	 "CGT":
				codonINFO[52] = codonINFO[52] + 1;
				break;
			case	 "CGC":
				codonINFO[53] = codonINFO[53] + 1;
				break;
			case	 "CGA":
				codonINFO[54] = codonINFO[54] + 1;
				break;
			case	 "CGG":
				codonINFO[55] = codonINFO[55] + 1;
				break;
			case	 "AGA":
				codonINFO[56] = codonINFO[56] + 1;
				break;
			case	 "AGG":
				codonINFO[57] = codonINFO[57] + 1;
				break;
			case	 "GGT":
				codonINFO[58] = codonINFO[58] + 1;
				break;
			case	 "GGC":
				codonINFO[59] = codonINFO[59] + 1;
				break;
			case	 "GGA":
				codonINFO[60] = codonINFO[60] + 1;
				break;
			case	 "GGG":
				codonINFO[61] = codonINFO[61] + 1;
				break;
			case	 "TAA":
				codonINFO[62] = codonINFO[62] + 1;
				break;
			case	 "TAG":
				codonINFO[63] = codonINFO[63] + 1;
				break;
			case	 "TGA":
				codonINFO[64] = codonINFO[64] + 1;
				break;
			default:
				codonINFO[0] = codonINFO[0] + 1;
				break;
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
		String info;
		info = "Index" + "\t" + "Headinfo" + "\t" + "Status" + "\t" + 
				"TTT_(F)\t" + 
				"TTC_(F)\t" + 
				"TTA_(L)\t" + 
				"TTG_(L)\t" + 
				"CTT_(L)\t" + 
				"CTC_(L)\t" + 
				"CTA_(L)\t" + 
				"CTG_(L)\t" + 
				"ATT_(I)\t" + 
				"ATC_(I)\t" + 
				"ATA_(I)\t" + 
				"ATG_(M)\t" + 
				"GTT_(V)\t" + 
				"GTC_(V)\t" + 
				"GTA_(V)\t" + 
				"GTG_(V)\t" + 
				"TCT_(S)\t" + 
				"TCC_(S)\t" + 
				"TCA_(S)\t" + 
				"TCG_(S)\t" + 
				"AGT_(S)\t" + 
				"AGC_(S)\t" + 
				"CCT_(P)\t" + 
				"CCC_(P)\t" + 
				"CCA_(P)\t" + 
				"CCG_(P)\t" + 
				"ACT_(T)\t" + 
				"ACC_(T)\t" + 
				"ACA_(T)\t" + 
				"ACG_(T)\t" + 
				"TAT_(Y)\t" + 
				"TAC_(Y)\t" + 
				"GCT_(A)\t" + 
				"GCC_(A)\t" + 
				"GCA_(A)\t" + 
				"GCG_(A)\t" + 
				"CAT_(H)\t" + 
				"CAC_(H)\t" + 
				"CAA_(Q)\t" + 
				"CAG_(Q)\t" + 
				"AAT_(N)\t" + 
				"AAC_(N)\t" + 
				"AAA_(K)\t" + 
				"AAG_(K)\t" + 
				"GAT_(D)\t" + 
				"GAC_(D)\t" + 
				"GAA_(E)\t" + 
				"GAG_(E)\t" + 
				"TGT_(C)\t" + 
				"TGC_(C)\t" + 
				"TGG_(W)\t" + 
				"CGT_(R)\t" + 
				"CGC_(R)\t" + 
				"CGA_(R)\t" + 
				"CGG_(R)\t" + 
				"AGA_(R)\t" + 
				"AGG_(R)\t" + 
				"GGT_(G)\t" + 
				"GGC_(G)\t" + 
				"GGA_(G)\t" + 
				"GGG_(G)\t" + 
				"TAA_(N/A)\t" + 
				"TAG_(N/A)\t" + 
				"TGA_(N/A)\n";
		write(outPath,info);
		readCDS(CDSPath);
	}
}
