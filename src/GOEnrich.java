import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.*;
import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

public class GOEnrich {
	public enum Type{ensembl,go};
	
	public static void main(String[] args) {
		
		String oboPath ="";
		String goname ="";
		String mappingPath ="";
		String simulationPath = "";
		
		String outputPath = "";
		String ooutputPath = "";
		
		int minsize = -1;
		int maxsize = -1;
		
		Type type = null;

		
		for(int i =0; i < args.length-1; i++){
			if(args[i].equals("-obo")){
				oboPath = args[i+1];
				i++;
			}
			else if(args[i].equals("-root")){
				goname = args[i+1];
				i++; 
			}
			else if(args[i].equals("-mapping")){
				mappingPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-mappingtype")){
				if("ensembl".equals(args[i+1])){
					type = Type.ensembl;
				}
				else if("go".equals(args[i+1])){
					type = Type.go;
				}
				i++; 
			}
			else if(args[i].equals("-enrich")){
				simulationPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-o")){
				outputPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-overlapout")){
				ooutputPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-minsize")){
				minsize = Integer.parseInt(args[i+1]);
				i++; 
			}
			else if(args[i].equals("-maxsize")){
				maxsize = Integer.parseInt(args[i+1]);
				i++; 
			}
		}
		
//		System.out.println(oboPath);
//		System.out.println(goname);
//		System.out.println(mappingPath);
//		System.out.println(simulationPath);
//		System.out.println(outputPath);
//		System.out.println(minsize);
//		System.out.println(maxsize);
		
		if(oboPath.equals("") || goname.equals("") || mappingPath.equals("") || simulationPath.equals("") || outputPath.equals("") || minsize == -1 || maxsize == -1){
			System.out.println("Usage Info:\n-obo <obo_file>\n-root <GO namespace>\n-mapping <gene2go_mapping>\n-mappingtype [ensembl|go]\n[-overlapout overlap_out_tsv]\n-enrich <simulation file>\n-o <output_tsv>\n-minsize <int>\n-maxsize <int>");
		}
		else{
			
			Path oboFilePath = Paths.get(oboPath);
			Path simulationFilePath = Paths.get(simulationPath);
			
			File oboFile = oboFilePath.toFile();			
			File simFile = simulationFilePath.toFile();
			try {
				GOEnrich goe = new GOEnrich(oboFile, goname, type, mappingPath, simFile, minsize, maxsize);
//				goe.getOutput(outputPath);
				if(!ooutputPath.equals("")){
					goe.getOverlapOut(ooutputPath);
				}
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	HashMap<String,GOEntry> go_entries;
	HashMap<String,HashSet<String>> overlaps;
	HashMap<String,Gene> genes;
	HashMap<String,HashSet<String>> gene_go;
	
	HashSet<String> enrichedGO;
	HashSet<String> enrichedGene;
	HashSet<String> signif;
	
	private int minSize;
	private int maxSize;
	
	
	public GOEnrich(File oboFile, String goname, Type mtype, String mappingPath, File simFile, int minSize, int maxSize) throws IOException{
		
		go_entries = new HashMap<String,GOEntry>();
		overlaps = new HashMap<String,HashSet<String>>();
		genes = new HashMap<String,Gene>();
		gene_go = new HashMap<String,HashSet<String>>();
		
		enrichedGO = new HashSet<String>();
		enrichedGene = new HashSet<String>();
		signif = new HashSet<String>();
		
		this.minSize = minSize;
		this.maxSize = maxSize;
		
		BufferedReader obobr = new BufferedReader (new FileReader(oboFile));
		String line ="";
		boolean inGOTerm = true;
		HashSet<String> goterm = new HashSet<String>();
    	String go_id = "";
    	String go_name = "";
    	String go_namespace ="";
    	HashSet<String> isa = null;
    	boolean obsolete = false;
	    while ((line = obobr.readLine()) != null){
	    	inGOTerm = !(line.contains("[T"));
	        if(inGOTerm){
	        	goterm.add(line);
	        }
	        else{
	        	isa = new HashSet<String>();
	        	for(String ti: goterm){
	        		if(ti.startsWith("id:")){
	        			String[] tisplit = ti.split(": ");
	        			go_id = tisplit[1];
//	        			if(go_id.equals("GO:1903034")){
//	        				System.out.println(goterm.toString());
//	        			}
	        		}
	        		else if(ti.startsWith("name:")){
	        			String[] tisplit = ti.split(": ");
	        			go_name = tisplit[1];
	        		}
	        		else if(ti.startsWith("namespace:")){
	        			String[] tisplit = ti.split(": ");
	        			go_namespace = tisplit[1];
	        		}
	        		else if(ti.startsWith("is_a:")){
	        			String[] tisplit = ti.split(": ");
	        			tisplit = tisplit[1].split(" ! ");
	        			isa.add(tisplit[0]);
	        		}
//	        		else if(ti.startsWith("is_obsolete:")){
//	        			String[] tisplit = ti.split(": ");
//	        			obsolete = Boolean.valueOf(tisplit[1]);
//	        		}
	        	}
	        	if(!go_id.equals("")){
//        			if(go_id.equals("GO:1903034")){
//        				System.out.println("I was here!");
//        			}
	        		if(go_namespace.equals(goname)){
		        		GOEntry newGOEntry = new GOEntry(go_id,go_name,isa);
		        		go_entries.put(go_id,newGOEntry);
//	        			obsolete = false;
	        		}
	        	}
	        	
	        	go_id = "";
	        	go_name = "";
	        	go_namespace ="";
	        	goterm.clear();
	        }
	    }
	    
//	    System.out.println(go_entries.containsKey("GO:1903034"));
	    
	    obobr.close();
	    
	    for(String go_entry: go_entries.keySet()){
	    	setGOEntry(go_entry);
	    }
	    line ="";
	    if(mtype==Type.ensembl){
			Path mappingFilePath = Paths.get(mappingPath);
			
			File mappingFile = mappingFilePath.toFile();
			BufferedReader mapbr = new BufferedReader (new FileReader(mappingFile));
			line = mapbr.readLine();
			while ((line = mapbr.readLine()) != null){
				String[] lineSplit = line.split("\t");
				String geneid = lineSplit[1];
				if(!geneid.equals("")){
					Gene curGene = new Gene(geneid);
					genes.put(geneid,curGene);
					if(!gene_go.containsKey(geneid)){
						gene_go.put(geneid, new HashSet<String>());
					}
					String[] goes = lineSplit[2].split("\\|");
					for(int i = 0; i < goes.length; i++){
//						System.out.println(goes[i]);
						if(go_entries.containsKey(goes[i])){
    						GOEntry goe = go_entries.get(goes[i]);
    						goe.addGene(geneid);
    						setDist(goes[i]);
    						for(String predgoeid : goe.getPred()){
    							go_entries.get(predgoeid).addGene(geneid);
//    							gene_go.get(geneid).add(predgoeid);
    						}
    						gene_go.get(geneid).add(goes[i]);
	    				}
					}
					
				}
			}
			mapbr.close();
	    }
	    else{
	    	FileInputStream gafinput = new FileInputStream(mappingPath);
	    	GZIPInputStream gafzip = new GZIPInputStream(gafinput);
	    	InputStreamReader ireader = new InputStreamReader(gafzip);
	    	BufferedReader gafbr = new BufferedReader(ireader);

	    	String geneid = "";
	    	String goeid ="";
	    	HashSet<String> curGOEntry = new HashSet<String>();
	    	while ((line = gafbr.readLine()) != null) {
	    		if(!line.startsWith("!")){
	    			String[] lineSplit = line.split("\t");
	    			if(lineSplit[3].equals("")){
	    				geneid = lineSplit[2];
	    				Gene newGene = new Gene(lineSplit[2]);
	    				genes.put(geneid,newGene);
						if(!gene_go.containsKey(geneid)){
							gene_go.put(geneid, new HashSet<String>());
						}
						goeid = lineSplit[4];
						if(go_entries.containsKey(goeid)){
    						GOEntry goe = go_entries.get(goeid);
    						goe.addGene(geneid);
    						setDist(goeid);
    						for(String predgoeid : goe.getPred()){
    							go_entries.get(predgoeid).addGene(geneid);
//    							gene_go.get(geneid).add(predgoeid);
    						}
    						gene_go.get(geneid).add(goeid);
	    				}
	    			}
	    		}
	    	}
			gafbr.close();
			ireader.close();
			gafzip.close();
			gafinput.close();
	    }
	    
	    for(String geneid : gene_go.keySet()){
	    	for(String goe1: gene_go.get(geneid)){
		    	for(String goe2: gene_go.get(geneid)){
		    		addOverlap(goe1,goe2);
		    	}
//		    	if(goe1.equals("GO:0015491")){
//		    		System.out.println(geneid);
//		    		System.out.println(gene_go.get(geneid).toString());
//		    	}
	    	}
	    }
	    
	    System.out.println(gene_go.get("SLC3A2").toString());
	    
	    System.out.println(go_entries.get("GO:0015171").toString());	    
	    System.out.println(go_entries.get("GO:0015491").toString());
	    
	    
//	    for( String go1 : overlaps.keySet()){
//	    	System.out.println(go1 + ": " + overlaps.get(go1).toString());
//	    }
	    
//	    for(String go_entry: go_entries.keySet()){
//	    	System.out.println(go_entries.get(go_entry).toString());
//	    }
	    
		BufferedReader simbr = new BufferedReader (new FileReader(simFile));
		line ="";

	    while ((line = simbr.readLine()) != null){
	    	if(line.startsWith("#")){
	    		enrichedGO.add(line.substring(1));
//	    		System.out.println(line.substring(1));
//	    		System.out.println(go_entries.containsKey(line.substring(1)));
//	    		System.out.println(go_entries.get(line.substring(1)).toString());
	    		
	    	}
	    	else if(line.equals("id\tfc\tsignif")){}
	    	else{
	    		String[] lineSplit = line.split("\t");
	    		String geneid = lineSplit[0];
	    		enrichedGene.add(geneid);
	    		if(genes.containsKey(geneid)){
	    			genes.get(geneid).setFC(Double.parseDouble(lineSplit[1]));
	    			if(lineSplit[2].equals("true")){
	    				signif.add(geneid);
	    			}
	    		}
	    	}
	    }
	    
	    simbr.close();
	}
	
	@SuppressWarnings("unchecked")
	public HashSet<String> setGOEntry(String id, HashSet<String> succs){
		GOEntry curGOE = go_entries.get(id);
		HashSet<String> parents = curGOE.is_a();
		succs.add(id);
		for(String parentId: parents){
			GOEntry parentGOE = go_entries.get(parentId);
			if(parentGOE  != null){
				parentGOE.addSucc(succs);
				HashSet<String> predIds = parentGOE.getPred();
				if(predIds.size() > 0){
					for(String predId : predIds){
						if(go_entries.get(predId)  != null){
							go_entries.get(predId).addSucc(succs);
						}
					}
					curGOE.addPred(predIds);
				}
				else{
					HashSet<String> newSuccs = (HashSet<String>) parentGOE.getSucc().clone();
					curGOE.addPred(setGOEntry(parentId, newSuccs));
				}
				curGOE.addPred(parentId);
			}
		};
		return curGOE.getPred();
	}
	
	@SuppressWarnings("unchecked")
	public void setGOEntry(String id){
		GOEntry curGOE = go_entries.get(id);
		if(curGOE.getPred().size() == 0){
		HashSet<String> parents = curGOE.is_a();
		HashSet<String> succs  = new HashSet<String>();
		succs.add(id);
		for(String parentId: parents){
			GOEntry parentGOE = go_entries.get(parentId);
			if(parentGOE  != null){
				parentGOE.addSucc(succs);
				HashSet<String> predIds = parentGOE.getPred();
				if(predIds.size() > 0){
					for(String predId : predIds){
						if(go_entries.get(predId)  != null){
							go_entries.get(predId).addSucc(succs);
						}
					}
					curGOE.addPred(predIds);
				}
				else{
					HashSet<String> newSuccs = new HashSet<String>(parentGOE.getSucc());
					curGOE.addPred(setGOEntry(parentId, newSuccs));
				}
				curGOE.addPred(parentId);
			}
		}
		}
	}
	
	public void addOverlap(String goe1, String goe2){
		boolean already_included = false;
		// goe2 before goe1
		if(goe1.compareTo(goe2) > 0){
			if(overlaps.containsKey(goe2)){
				already_included = overlaps.get(goe2).contains(goe1);
			}
		}
		else{
			if(overlaps.containsKey(goe1)){
				already_included = overlaps.get(goe1).contains(goe2);
			}
		}
		if(!already_included){
			HashSet<String> set1 = new HashSet<String>(go_entries.get(goe1).getPred());
			set1.add(goe1);
			HashSet<String> set2 = new HashSet<String>(go_entries.get(goe2).getPred());
			set2.add(goe2);
			
			for(String ov1: set1){
				for(String ov2: set2){
					if(ov1.compareTo(ov2) > 0){
						if(!overlaps.containsKey(ov2)){
							overlaps.put(ov2,new HashSet<String>());
						}
						overlaps.get(ov2).add(ov1);
					}
					else{
						if(!overlaps.containsKey(ov1)){
							overlaps.put(ov1,new HashSet<String>());
						}
						overlaps.get(ov1).add(ov2);						
					}
				}
			}
		}
		
	}
	
	public HashMap<String,Integer> setDist(String goeid){
		
		GOEntry goe = go_entries.get(goeid);
		HashMap<String,Integer> stor = new HashMap<String,Integer>();
		
		
		if(goe.getDistance().isEmpty()){
			for(String parent_id : goe.is_a()){
				if(go_entries.containsKey(parent_id)){
					stor = this.setDist(parent_id);
				}
				else{
					stor = new HashMap<String,Integer>();
				}
				goe.addDistance(stor, parent_id);
			}
			goe.putDistance(goeid,0);
		}
		return goe.getDistance();
	}
	
	public void getOverlapOut(String outputPath) throws IOException {
		String tab = "\t";
		String brk = "\n";
		char dop = ':';
		char sip = '|';
		char lin = '-';
		char com = ',';
		char spa = ' ';
		
		StringBuilder resultBuilder = new StringBuilder("");
		resultBuilder.append("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent\n");
		FileWriter outputWriter = new FileWriter(outputPath,false);
		outputWriter.write(resultBuilder.toString());
		resultBuilder.setLength(0);
		
		ArrayList<String> overlapKeyList = new ArrayList<String>(overlaps.keySet());
		overlapKeyList.sort(new StringComparator());
		
		GOEntry goe1;
		GOEntry goe2;
		
		HashSet<String> goe1Pred;
		HashSet<String> goe1Succ;
		HashSet<String> goe1Gene;
		HashMap<String,Integer> goe1Dist; 
		int goe1Size;
		
		HashSet<String> goe2Pred;
//		HashSet<String> goe2Succ;
		HashSet<String> goe2Gene;
		HashMap<String,Integer> goe2Dist;
		int goe2Size;
		
		boolean relative = false;
		int dist = -1;
		int overlap = -1;
		double percent = -1;
		

		
		for(String goe1Id: overlapKeyList){
			goe1 = go_entries.get(goe1Id);
			ArrayList<String> goe2List = new ArrayList<String>(overlaps.get(goe1Id));
			goe2List.sort(new StringComparator());
			for(String goe2Id: goe2List){
				if(!goe2Id.equals(goe1Id)){
					goe2 = go_entries.get(goe2Id);
					
					goe1Pred = new HashSet<String>(goe1.getPred());
					goe1Succ = new HashSet<String>(goe1.getSucc());
					goe1Gene = new HashSet<String>(goe1.getGenes());
					goe1Dist = new HashMap<String,Integer>(goe1.getDistance());
					goe1Size = goe1Gene.size();
					
					goe2Pred = new HashSet<String>(goe2.getPred());
//					goe2Succ = new HashSet<String>(goe2.getSucc());
					goe2Gene = new HashSet<String>(goe2.getGenes());
					goe2Dist = new HashMap<String,Integer>(goe2.getDistance());
					goe2Size = goe2Gene.size();
					
					if(goe1Id.equals("GO:0015171") && goe2Id.equals("GO:0015491")){
						System.out.println(goe1Size);
						System.out.println(goe2Size);
					}
				    
					
					if(goe1Size >= minSize && goe2Size >= minSize && goe1Size <= maxSize && goe2Size <= maxSize){
						resultBuilder.append(goe1Id);
						resultBuilder.append(tab);
						resultBuilder.append(goe2Id);
						resultBuilder.append(tab);
						
						relative = goe1Pred.contains(goe2Id) || goe1Succ.contains(goe2Id);
						
						resultBuilder.append(relative);
						resultBuilder.append(tab);
						
						if(relative){
							if(goe1Dist.containsKey(goe2Id)){
								dist = goe1Dist.get(goe2Id);
							}
							else if(goe2Dist.containsKey(goe1Id)){
								dist = goe2Dist.get(goe1Id);
							}
						}
						else{
							goe1Pred.retainAll(goe2Pred);
							for(String compred : goe1Pred){
								if(dist == -1){
									dist = goe1Dist.get(compred) + goe2Dist.get(compred);
								}
								else{
									dist = Math.min(dist, goe1Dist.get(compred) + goe2Dist.get(compred));
								}
							}
						}
						resultBuilder.append(dist);
						resultBuilder.append(tab);
						
						goe1Gene.retainAll(goe2Gene);
						overlap = goe1Gene.size();
						
						resultBuilder.append(overlap);
						resultBuilder.append(tab);
						
						percent = (overlap* 100.0)/goe1Size;
						percent = Math.max(percent, (overlap* 100.0)/goe2Size);
						
						resultBuilder.append(percent);
						resultBuilder.append(brk);
						
						outputWriter.write(resultBuilder.toString());
						resultBuilder.setLength(0);
					}
				}
			}
		}
		outputWriter.close();
		
	}

	public void getOutput(String outputPath) throws IOException {
		String tab = "\t";
		String brk = "\n";
		char dop = ':';
		char sip = '|';
		char lin = '-';
		char com = ',';
		char spa = ' ';
		char str = '*';

		StringBuilder resultBuilder = new StringBuilder("");
		resultBuilder.append("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true\n");
		FileWriter outputWriter = new FileWriter(outputPath,false);
		outputWriter.write(resultBuilder.toString());
		resultBuilder.setLength(0);
		
		ArrayList<GOEntry> outputGOE = new ArrayList<GOEntry>();
		ArrayList<Double> hgValues = new ArrayList<Double>();
		ArrayList<Double> fejValues = new ArrayList<Double>();
		ArrayList<Double> kssValues = new ArrayList<Double>();
		ArrayList<Double> kspValues = new ArrayList<Double>();
		
		GOEntry goe;
		HashSet<String> goeGene;
		int size;
		int noverlap;
		
		HypergeometricDistribution hd;
		HypergeometricDistribution fej;
		KolmogorovSmirnovTest kst;
		
		ArrayList<Double> backgroundFC = new ArrayList<Double>();
		
		enrichedGene.retainAll(genes.keySet());
		signif.retainAll(genes.keySet());
		int eSize = enrichedGene.size();
		int sSize = signif.size();
		
		enrichedGene.retainAll(genes.keySet());
		for(String eg : enrichedGene){
			backgroundFC.add(genes.get(eg).getFC());
		}
		
		double[] bfc  = new double[backgroundFC.size()];
		for(int i = 0; i < bfc.length; i++){
			bfc[i] = backgroundFC.get(i);
		}
		double[] fc;

		
		for(String goeid: go_entries.keySet()){
			goe = go_entries.get(goeid);
			goeGene = new HashSet<String>(goe.getGenes());
//			if(goeid.equals("GO:0006200")){
//				System.out.println(goe.toString());
//			}
			if(goeGene.size() >= minSize && goeGene.size() <= maxSize){
				goeGene.retainAll(enrichedGene);
				size = goeGene.size();
				goe.setSize(size);
				
				fc = new double[goeGene.size()];
				int j = 0;
				for(String g: goeGene){
					fc[j] = genes.get(g).getFC();
					j++;
				}
				
				goeGene.retainAll(signif);
				noverlap = goeGene.size();
				goe.setNOverlap(noverlap);
				
				goe.setTruth(enrichedGO.contains(goeid));
				outputGOE.add(goe);
				
				hd = new HypergeometricDistribution(eSize, sSize, size);
				
				hgValues.add(hd.upperCumulativeProbability(noverlap));
				
				fej = new HypergeometricDistribution(eSize-1, sSize-1, size-1);
				fejValues.add(fej.upperCumulativeProbability(noverlap-1));
				
				kst = new KolmogorovSmirnovTest();
				kssValues.add(kst.kolmogorovSmirnovStatistic(fc, bfc));
				kspValues.add(kst.kolmogorovSmirnovTest(fc, bfc));
			}
		}
		
		Double[] hgArray = new Double[hgValues.size()];
		Double[] fejArray = new Double[fejValues.size()];
		Double[] kspArray = new Double[kspValues.size()];
		
		double[] hghbValues = BenjaminiHochberg(hgValues.toArray(hgArray));
		double[] fejhbValues = BenjaminiHochberg(fejValues.toArray(fejArray));
		double[] kshbValues = BenjaminiHochberg(kspValues.toArray(kspArray));
		
		StringBuilder spBuilder = new StringBuilder();
		
		int minDist = -1;
		int dist;
		GOEntry minGOE = null;
		GOEntry trueGOE = null;
		
		HashSet<String> goePred = null;
		HashMap<String, Integer> goe1Dist;
		HashMap<String, Integer> goe2Dist;
		HashSet<GOEntry> trueGO = new HashSet<GOEntry>();
		for(String eGOE:enrichedGO){
			trueGO.add(go_entries.get(eGOE));
//			System.out.println(go_entries.get(eGOE).getName());
		}
		
		for(int i = 0; i < outputGOE.size(); i++){
			goe = outputGOE.get(i);
			resultBuilder.append(goe.getId());
			resultBuilder.append(tab);
			
			resultBuilder.append(goe.getName());
			resultBuilder.append(tab);
			
			resultBuilder.append(goe.getSize());
			resultBuilder.append(tab);
			
			resultBuilder.append(goe.getTruth());
			resultBuilder.append(tab);
			
			resultBuilder.append(goe.getNOverlap());
			resultBuilder.append(tab);
			
			resultBuilder.append(hgValues.get(i));
			resultBuilder.append(tab);
			resultBuilder.append(hghbValues[i]);
			resultBuilder.append(tab);
			
			resultBuilder.append(fejValues.get(i));
			resultBuilder.append(tab);
			resultBuilder.append(fejhbValues[i]);
			resultBuilder.append(tab);
			
			resultBuilder.append(kssValues.get(i));
			resultBuilder.append(tab);
			resultBuilder.append(kspValues.get(i));
			resultBuilder.append(tab);
			resultBuilder.append(kshbValues[i]);
			
			if(!goe.getTruth()){
				goe1Dist = goe.getDistance();
				if(goe1Dist.isEmpty()){
					setDist(goe.getId());
				}
				minDist = -1;
				minGOE = null;
				trueGOE = null;
				for(GOEntry tGOE : trueGO){
					goe2Dist = tGOE.getDistance();
					if(goe1Dist.containsKey(tGOE.getId())){
						if(minDist == -1){
							minDist = goe1Dist.get(tGOE.getId());
							minGOE = tGOE;
							trueGOE = tGOE;
						}
						else{
							dist = goe1Dist.get(tGOE.getId());
							if(minDist > dist){
								minDist = dist;
								minGOE = tGOE;
								trueGOE = tGOE;						
							}
						}
					}
					else if(goe2Dist.containsKey(goe.getId())){
						if(minDist == -1){
							minDist = goe2Dist.get(goe.getId());
							minGOE = goe;
							trueGOE = tGOE;
						}
						else{
							dist = goe2Dist.get(goe.getId());
							if(minDist > dist){
								minDist = dist;
								minGOE = goe;
								trueGOE = tGOE;						
							}
						}
					}
					else{
						goePred = new HashSet<String>(goe.getPred());
						goePred.retainAll(tGOE.getPred());
						for(String compred : goePred){
							if(minDist == -1){
								minDist = goe1Dist.get(compred) + goe2Dist.get(compred);
								minGOE = go_entries.get(compred);
								trueGOE = tGOE;
							}
							else{
								dist = goe1Dist.get(compred) + goe2Dist.get(compred);
								if(minDist > dist){
									minDist = dist;
									minGOE = go_entries.get(compred);
									trueGOE = tGOE;									
								}
							}
						}
					}
				}
				
//				System.out.println(goe.toString());
				
				if(!(minGOE != null)){
				}else{
				
				ArrayList<String> goeToLCA =  getPath(goe,minGOE,goe.getDistance().get(minGOE.getId()));
				
				spBuilder.append(go_entries.get(goeToLCA.get(0)).getName());
				for(int j = 1; j < goeToLCA.size(); j++){
					spBuilder.append(sip);
					spBuilder.append(go_entries.get(goeToLCA.get(j)).getName());
				}
				spBuilder.append(str);
				
				ArrayList<String> tgoeToLCA =  getPath(trueGOE,minGOE,trueGOE.getDistance().get(minGOE.getId()));
				for(int j = tgoeToLCA.size()-2; j > 0; j--){
					spBuilder.append(sip);
					spBuilder.append(go_entries.get(tgoeToLCA.get(j)).getName());					
				}
				}
				resultBuilder.append(tab);
				resultBuilder.append(spBuilder.toString());
				spBuilder.setLength(0);
			}
			
			
			resultBuilder.append(brk);
			
			outputWriter.write(resultBuilder.toString());
			resultBuilder.setLength(0);
		}
		
	}
	
	public double[] BenjaminiHochberg(Double[] pValues) {
		
		int length = pValues.length;
        Arrays.sort(pValues);
        double[] apValues = new double[length];

        for (int i = length - 1; i >= 0; i--) {
            if (i == length - 1) {
                apValues[i] = pValues[i];
            } else {
                double l = apValues[i + 1];
                double r = (length / (double) (i+1)) * pValues[i];
                apValues[i] = Math.min(l, r);
            }
        }
        return apValues;
	}
	
	public ArrayList<String> getPath(GOEntry goeFrom, GOEntry goeTo, int length){
		ArrayList<String> result = null;
		String goeToId = goeTo.getId();
		String goeFromId = goeFrom.getId();
		GOEntry parent;
		if(goeFrom.getPath(goeToId) != null){
			return goeFrom.getPath(goeToId);
		}
		else if(goeFromId.equals(goeToId)){
			result = new ArrayList<String>();
			result.add(goeFrom.getId());
			goeFrom.setPath(result);
		}
		else{
			for(String goePId : goeFrom.is_a()){
				parent = go_entries.get(goePId);
				if(parent.getDistance().containsKey(goeToId)){
					if(parent.getDistance().get(goeToId) == length-1){
						if(parent.getPath(goeToId) == null){
							result = getPath(parent,goeTo,length-1);
						}
						else{
							result = parent.getPath(goeToId);
						}
						result.add(0, goeFromId);
						goeFrom.setPath(result);
						break;
					}
				}
			}
		}
		return result;
	}
	
	class StringComparator implements Comparator<String>
	{
	    public int compare(String x1, String x2)
	    {
	        return x1.compareTo(x2);
	    }
	}
	
}
