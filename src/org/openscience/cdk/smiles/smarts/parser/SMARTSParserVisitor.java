package org.openscience.cdk.smiles.smarts.parser;


//FIXMED CLASE VACIA
public interface SMARTSParserVisitor {

	  public Object visit(SimpleNode node, Object data);
	  public Object visit(ASTStart node, Object data);
	  public Object visit(ASTReaction node, Object data);
	  public Object visit(ASTGroup node, Object data);
	  public Object visit(ASTSmarts node, Object data);
	  public Object visit(ASTAtom node, Object data);
	  public Object visit(ASTLowAndBond node, Object data);
	  public Object visit(ASTOrBond node, Object data);
	  public Object visit(ASTExplicitHighAndBond node, Object data);
	  public Object visit(ASTImplicitHighAndBond node, Object data);
	  public Object visit(ASTNotBond node, Object data);
	  public Object visit(ASTSimpleBond node, Object data);
	  public Object visit(ASTExplicitAtom node, Object data);
	  public Object visit(ASTLowAndExpression node, Object data);
	  public Object visit(ASTOrExpression node, Object data);
	  public Object visit(ASTExplicitHighAndExpression node, Object data);
	  public Object visit(ASTImplicitHighAndExpression node, Object data);
	  public Object visit(ASTNotExpression node, Object data);
	  public Object visit(ASTRecursiveSmartsExpression node, Object data);
	  public Object visit(ASTTotalHCount node, Object data);
	  public Object visit(ASTImplicitHCount node, Object data);
	  public Object visit(ASTExplicitConnectivity node, Object data);
	  public Object visit(ASTAtomicNumber node, Object data);
	  public Object visit(ASTHybrdizationNumber node, Object data);
	  public Object visit(ASTCharge node, Object data);
	  public Object visit(ASTRingConnectivity node, Object data);
	  public Object visit(ASTPeriodicGroupNumber node, Object data);
	  public Object visit(ASTTotalConnectivity node, Object data);
	  public Object visit(ASTValence node, Object data);
	  public Object visit(ASTRingMembership node, Object data);
	  public Object visit(ASTSmallestRingSize node, Object data);
	  public Object visit(ASTAliphatic node, Object data);
	  public Object visit(ASTNonCHHeavyAtom node, Object data);
	  public Object visit(ASTAromatic node, Object data);
	  public Object visit(ASTAnyAtom node, Object data);
	  public Object visit(ASTAtomicMass node, Object data);
	  public Object visit(ASTRingIdentifier node, Object data);
	  public Object visit(ASTChirality node, Object data);
	  public Object visit(ASTElement node, Object data);
	

}
