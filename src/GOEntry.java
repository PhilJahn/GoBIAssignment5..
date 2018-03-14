import java.util.HashSet;

public class GOEntry {
	
	private String id;
	private String name;
	private HashSet<String> isa;
	private HashSet<String> pred;
	private HashSet<String> succ;
	private HashSet<String> gene;
	
	public GOEntry(String id, String name, HashSet<String> isa){
		this.id = id;
		this.name = name;
		this.isa = isa;
	
		pred = new HashSet<String>();
		succ = new HashSet<String>();
		gene = new HashSet<String>();
	}
	
	public String getId(){
		return id;
	}
	
	public String getName(){
		return name;
	}
	
	public HashSet<String> is_a(){
		return isa;
	}
	
	public HashSet<String> getPred(){
		return pred;
	}
	
	public HashSet<String> getSucc(){
		return succ;
	}
	
	public HashSet<String> getGenes(){
		return gene;
	}
	
	public void addGene(String gene){
		this.gene.add(gene);
	}
	
	public void addSucc(HashSet<String> succs){
		succ.addAll(succs);
	}
	
	public void addSucc(String succ){
		this.succ.add(succ);
	}
	
	public void addPred(HashSet<String> preds){
		pred.addAll(preds);
	}
	
	public void addPred(String pred){
		this.pred.add(pred);
	}
	
	public String toString(){
		String brk ="\n";
		String tab = "\t";
		StringBuilder entryBuilder = new StringBuilder("Id:");
		entryBuilder.append(tab);
		entryBuilder.append(id);
		entryBuilder.append(brk);
		entryBuilder.append("Name:");
		entryBuilder.append(tab);
		entryBuilder.append(name);
		entryBuilder.append(brk);
		entryBuilder.append("is a:");
		entryBuilder.append(tab);
		entryBuilder.append(isa.toString());
		entryBuilder.append(brk);
		entryBuilder.append("Pred:");
		entryBuilder.append(tab);
		entryBuilder.append(pred.toString());
		entryBuilder.append(brk);
		entryBuilder.append("Succ:");
		entryBuilder.append(tab);
		entryBuilder.append(succ.toString());
		entryBuilder.append(brk);
		entryBuilder.append("Genes:");
		entryBuilder.append(tab);
		entryBuilder.append(gene.toString());
		
		return entryBuilder.toString();
	}

}
