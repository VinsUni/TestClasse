/*
 * $Id: DistanceAlignment.java 1630 2011-09-15 20:29:40Z euzenat $
 *
 * Copyright (C) INRIA, 2003-2011
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 */

package fr.inrialpes.exmo.align.impl; 

import java.lang.ClassCastException;
import java.util.TreeSet;
import java.util.SortedSet;
import java.util.Comparator;
import java.util.Properties;

import org.semanticweb.owl.align.Alignment;
import org.semanticweb.owl.align.AlignmentProcess;
import org.semanticweb.owl.align.AlignmentException;
import org.semanticweb.owl.align.Cell;

import fr.inrialpes.exmo.ontowrap.OntowrapException;

import fr.inrialpes.exmo.ontosim.util.HungarianAlgorithm;

/**
 * The mother class for distance or similarity-based alignments.
 * It is abstract because it does not provide an implemented similarity measure
 * Otherwise everything is fine.
 *
 * This class should work with similarity and distances, as soon as, the used
 * similarity structure is defined as such.
 *
 * @author Jérôme Euzenat
 * @version $Id: DistanceAlignment.java 1630 2011-09-15 20:29:40Z euzenat $ 
 */

public abstract class DistanceAlignment extends fr.inrialpes.exmo.align.impl.ObjectAlignment implements AlignmentProcess {
    Similarity sim;

    /** Creation **/
    public DistanceAlignment() {};

    public void setSimilarity( Similarity s ) { sim = s; }
    public Similarity getSimilarity() { return sim; }

    /**
     * Process matching
     * - create distance data structures,
     * - compute distance or similarity
     * - extract alignment
     **/
    public void align( Alignment alignment, Properties params ) throws AlignmentException {
	loadInit( alignment );
	if (  params.getProperty("type") != null ) 
	    setType( params.getProperty("type") );
	// This is a 1:1 alignment in fact
	else setType("11");
	if ( sim == null )
	    throw new AlignmentException("DistanceAlignment: requires a similarity measure");

	sim.initialize( ontology1(), ontology2(), init );
	sim.compute( params );
	if ( params.getProperty("printMatrix") != null ) printDistanceMatrix(params);
	extract( getType(), params );
    }

    /**
     * Prints the distance matrix
     */
    public void printDistanceMatrix( Properties params ){
	String algName = params.getProperty("algName");
	String metric = "distance";
	if ( sim.getSimilarity() ) metric = "similarity";
	if ( algName == null ) algName = getClass().toString();
	System.out.println("\\documentclass{article}\n");
	System.out.println("\\usepackage{graphics}\n");
	System.out.println("\\begin{document}\n");
	System.out.println("\\begin{table}");
	sim.printClassSimilarityMatrix("tex");
	System.out.println("\\caption{Class "+metric+" with measure "+algName+".}" );
	System.out.println("\\end{table}");
	System.out.println();
	System.out.println("\\begin{table}");
	sim.printPropertySimilarityMatrix("tex");
	System.out.println("\\caption{Property "+metric+" with measure "+algName+".}" );
	System.out.println("\\end{table}");
	System.out.println();
	System.out.println("\\begin{table}");
	sim.printIndividualSimilarityMatrix("tex");
	System.out.println("\\caption{Individual "+metric+" with measure "+algName+".}" );
	System.out.println("\\end{table}");
	System.out.println("\n\\end{document}");
    }

    /**
     * Suppresses the distance matrix
     */
    public void cleanUp() {
	sim = null;
    }

    /**
     * Extract the alignment form the Similarity
     * There are theoretically 16 types of extractors composing the
     * characteristics
     * [q]estion mark = ?, one or zero relation
     * [s]tar = *, one, zero or many relations
     * [1] = 1, exactly one relation
     * [p]lus = +, one or many relations
     * for each place of the relation. Howerver, since it is not possible from a matrics to guarantee that one object will be in at least one relation, this is restricted to the four following types:
     * ?? (covering 11, 1? and ?1)
     * ** (covering ++, *+ and +*)
     * ?* (covering 1*, 1+ and ?+)
     * *? (covering +?, *1 and +1)
     */

    public Alignment extract(String type, Properties params) throws AlignmentException {
	double threshold = 0.;
	if (  params.getProperty("threshold") != null )
	    threshold = Double.parseDouble( params.getProperty("threshold") );

	//System.err.println("The type is "+type+" with length = "+type.length());
	return this.extractMethodA(type, params, threshold);
    }

	/**
	 * New Method extractMethodA()
	 */
	private Alignment extractMethodA(String type, Properties params, double threshold) throws AlignmentException {
		if ( type.equals("?*") || type.equals("1*") || type.equals("?+") || type.equals("1+") ) {
			return extractqs(threshold, params);
		} else{
			return this.extractMethodA1(type, params, threshold);
		}
	}
	// JE: It is now certainly possible to virtualise extraction as it has
    // been done for printing matrix in MatrixMeasure (todo)

	/**
	 * New Method extractMethodA1()
	 */
	private Alignment extractMethodA1(String type, Properties params, double threshold) throws AlignmentException {
		if ( type.equals("??") || type.equals("1?") || type.equals("?1") || type.equals("11") ) {
			return extractqq(threshold, params);
		}else{
			return extractMethodA2(type, params, threshold);
		}
	}

	/**
	 * New Method extractMethodA2()
	 */
	private Alignment extractMethodA2(String type, Properties params, double threshold) throws AlignmentException {
		if ( type.equals("*?") || type.equals("+?") || type.equals("*1") || type.equals("+1") ) {
			return extractqs(threshold, params);
		}else{
			return extractMethodA3(type, params, threshold);
		}
	}

	/**
	 * New Method extractMethodA3()
	 */
	private Alignment extractMethodA3(String type, Properties params, double threshold ) throws AlignmentException {
		if ( type.equals("**") || type.equals("+*") || type.equals("*+") || type.equals("++") ) {
			return extractss(threshold, params);
			// The else should be an error message
		}else{
			return extractMethodA4(type, params, threshold);
		}
	}

	/**
	 * New Method extractMethodA4()
	 */
	private Alignment extractMethodA4(String type, Properties params, double threshold) throws AlignmentException {
		if ( type.equals("**") || type.equals("+*") || type.equals("*+") || type.equals("++") ) {
			return extractss(threshold, params);
			// The else should be an error message
		}else {
			throw new AlignmentException("Unknown alignment type: " + type);
		}
	}


    /**
     * Extract the alignment of a ?* type
     * Non symmetric: for each entity of onto1, take the highest if superior to threshold
     * Complexity: O(n^2)
     */
    @SuppressWarnings("unchecked") //ConcatenatedIterator
    public Alignment extractqs( double threshold, Properties params) {
      double max = 0.;
      boolean found = false;
      double val = 0.;

      try {
	  // Extract for properties
	  ConcatenatedIterator pit1 = new
	      ConcatenatedIterator(ontology1().getObjectProperties().iterator(),
				   ontology1().getDataProperties().iterator());

	  this.extractqsForMethodA(pit1, found, max, threshold, val);

	  // Extract for classes
	  this.extractqsForMethodB(found, max, threshold, val);

	  // Extract for individuals
	  this.extractqsIfElseMethod(params, found, max, val, threshold);

      } catch (OntowrapException owex) { owex.printStackTrace(); //}
      } catch (AlignmentException alex) { alex.printStackTrace(); }
      return((Alignment)this);
    }

	/**
	 * New Method extractqsForMethodA()
	 */

	private void extractqsForMethodA(ConcatenatedIterator pit1, boolean found, double max, double threshold, double val ){
		for( Object prop1 : pit1 ){
			found = false; max = threshold; val = 0.;
			Object prop2 = null;
			ConcatenatedIterator pit2 = new
					ConcatenatedIterator(ontology2().getObjectProperties().iterator(),
					ontology2().getDataProperties().iterator());
			this.extractqsInnerForMethodA(pit2, double val, double max, Object prop1, Object prop2, boolean found);
			if ( found ) addAlignCell(prop1,prop2, "=", max);
		}
	}

	/**
	 * New Method extractqsInnerForMethodA()
	 */
	private void extractqsInnerForMethodA(ConcatenatedIterator pit2, double val, double max, Object prop1, Object prop2, boolean found){
		for ( Object current : pit2 ){
			if ( sim.getSimilarity() ) val = sim.getPropertySimilarity(prop1,current);
			else val =  1. - sim.getPropertySimilarity(prop1,current);
			if ( val > max) {
				found = true; max = val; prop2 = current;
			}
		}
	}

	/**
	 * New Method extractqsForMethodB()
	 */
	private void extractqsForMethodB(boolean found, double max, double threshold, double val){
		for ( Object class1 : ontology1().getClasses() ) {
			found = false; max = threshold; val = 0;
			Object class2 = null;
			this.extractqsInnerForMethodB(class1, class2, found, max, val);
			if ( found ) addAlignCell(class1, class2, "=", max);
		}
	}

	/**
	 * New Method extractqsInnerForMethodB()
	 */
	private void extractqsInnerForMethodB(Object class1, Object class2, boolean found, double max, double val){
		for ( Object current : ontology2().getClasses() ) {
			if ( sim.getSimilarity() ) val = sim.getClassSimilarity(class1,current);
			else val = 1. - sim.getClassSimilarity(class1,current);
			if (val > max) {
				found = true; max = val; class2 = current;
			}
		}
	}

	/**
	 * New Method extractqsIfElseMethod()
	 */
	private void extractqsIfElseMethod(Properties params, boolean found, double max, double val, double threshold){
		if (  params.getProperty("noinst") == null ){
			for ( Object ind1 : ontology1().getIndividuals() ) {
				this.extractqsIfElseInnerMethodA(ind1, found, max, threshold, val );
			}
		}
	}

	/**
	 * New Method extractqsIfElseInnerMethodA()
	 */
	private void extractqsIfElseInnerMethodA(Object ind1, boolean found, double max, double val, double threshold){
		if ( ontology1().getEntityURI( ind1 ) != null ) {
			found = false; max = threshold; val = 0;
			Object ind2 = null;
			for ( Object current : ontology2().getIndividuals() ) {
				this.extractqsIfElseInnerMethodB(ind1, current, val, max, found);
			}
			if ( found ) addAlignCell(ind1,ind2, "=", max);
		}
	}

	/**
	 * New Method extractqsIfElseInnerMethodB()
	 */
	private void extractqsIfElseInnerMethodB(Object ind1, Object current, double val, double max, boolean found ){
		if ( ontology2().getEntityURI( current ) != null ) {
			if ( sim.getSimilarity() ) val = sim.getIndividualSimilarity( ind1, current );
			else val = 1 - sim.getIndividualSimilarity( ind1, current );
			if (val > max) {
				found = true; max = val; ind2 = current;
			}
		}
	}

	/**
     * Extract the alignment of a ** type
     * Symmetric: return all elements above threshold
     * Complexity: O(n^2)
     */
    @SuppressWarnings("unchecked") //ConcatenatedIterator
    public Alignment extractss( double threshold, Properties params) {
	double val = 0.;
	try {
	    // Extract for properties
	    ConcatenatedIterator pit1 = new 
		ConcatenatedIterator(ontology1().getObjectProperties().iterator(),
				     ontology1().getDataProperties().iterator());
		this.extractssForMethodA(pit1, val, threshold);

	    // Extract for classes
		this.extractssForMethodB(val, threshold);

	    // Extract for individuals
	    if (  params.getProperty("noinst") == null ){
			this.extractssForMethodC(val, threshold);

	    }
	} catch (OntowrapException owex) { owex.printStackTrace(); //}
	} catch (AlignmentException alex) { alex.printStackTrace(); }
	return((Alignment)this);
    }

	/**
	 * New Method extractssForMethodA()
	 */
	private void extractssForMethodA(ConcatenatedIterator pit1, double val, double threshold){
		for( Object prop1 : pit1 ){
			ConcatenatedIterator pit2 = new
					ConcatenatedIterator(ontology2().getObjectProperties().iterator(),
					ontology2().getDataProperties().iterator());
			this.extractssInnerForMethodA();

		}
	}

	/**
	 * New Method extractssInnerForMethodA()
	 */
	private void extractssInnerForMethodA(ConcatenatedIterator pit2, double val, double threshold){
		for ( Object prop2 : pit2 ){
			if ( sim.getSimilarity() ) val = sim.getPropertySimilarity(prop1,prop2);
			else val =  1. - sim.getPropertySimilarity(prop1,prop2);
			if ( val > threshold ) addAlignCell(prop1,prop2, "=", val);
		}
	}


	/**
	 * New Method extractssForMethodB()
	 */
	private void extractssForMethodB(double val, double threshold){
		for ( Object class1 : ontology1().getClasses() ) {
			for ( Object class2 : ontology2().getClasses() ) {
				if ( sim.getSimilarity() ) val = sim.getClassSimilarity(class1,class2);
				else val = 1. - sim.getClassSimilarity(class1,class2);
				if (val > threshold ) addAlignCell(class1, class2, "=", val);
			}
		}
	}

	/**
	 * New Method extractssForMethodC()
	 */
	private void extractssForMethodC(double val, double threshold){
		for ( Object ind1 : ontology1().getIndividuals() ) {
			if ( ontology1().getEntityURI( ind1 ) != null ) {
				this.extractssInnerForMethodC(val, threshold);
			}
		}
	}

	/**
	 * New Method extractssInnerForMethodC()
	 */
	private void extractssInnerForMethodC(double val, double threshold){
		for ( Object ind2 : ontology2().getIndividuals() ) {
			if ( ontology2().getEntityURI( ind2 ) != null ) {
				if ( sim.getSimilarity() ) val = sim.getIndividualSimilarity( ind1, ind2 );
				else val = 1 - sim.getIndividualSimilarity( ind1, ind2 );
				if ( val > threshold ) addAlignCell(ind1,ind2, "=", val);
			}
		}
	}


	/**
     * Extract the alignment of a ?? type
     * 
     * exact algorithm using the Hungarian method.
     * This algorithm contains several guards to prevent the HungarianAlgorithm to
     * raise problems:
     * - It invert column and rows when nbrows > nbcol (Hungarian loops)
     * - It prevents to generate alignments when one category has no elements.
     */
    @SuppressWarnings("unchecked") //ConcatenatedIterator
    public Alignment extractqq( double threshold, Properties params) {
	try {
	    // A STRAIGHTFORWARD IMPLEMENTATION
	    // (redoing the matrix instead of getting it)
	    // For each kind of stuff (cl, pr, ind)
	    // Create a matrix
	    int nbclasses1 = ontology1().nbClasses();
	    int nbclasses2 = ontology2().nbClasses();
	   this.extractqqIfElseMethodA(nbclasses1, nbclasses2, threshold);

	} catch ( AlignmentException alex) { alex.printStackTrace(); }
	catch ( OntowrapException owex) { owex.printStackTrace(); }
	// For properties
	try{
	    int nbprop1 = ontology1().nbProperties();
	    int nbprop2 = ontology2().nbProperties();
	   this.extractqqIfElseMethodB( nbprop1, nbprop2, threshold);

	} catch (AlignmentException alex) { alex.printStackTrace(); }
	catch (OntowrapException owex) { owex.printStackTrace(); }
	// For individuals
		this.extractqqIfElseMethodC(params, threshold);

	return((Alignment)this);
    }

    public int[][] callHungarianMethod( double[][] matrix, int i, int j ) {
	boolean transposed = false;
	if ( i > j ) { // transposed aray (because rows>columns).
	    matrix = HungarianAlgorithm.transpose(matrix);
	    transposed = true;
	}
	int[][] result = HungarianAlgorithm.hgAlgorithm( matrix, "max" );
	if ( transposed ) {
	    for( int k=0; k < result.length ; k++ ) { 
		int val = result[k][0]; result[k][0] = result[k][1]; result[k][1] = val; 
	    }
	    
	}
	return result;
    }

	/**
	 *
	 * @param nbclasses1
	 * @param nbclasses2
	 * @param threshold
	 */
    private void extractqqIfElseMethodA(int nbclasses1, int nbclasses2, double threshold){
		if ( nbclasses1 != 0 && nbclasses2 != 0 ) {
			double[][] matrix = new double[nbclasses1][nbclasses2];
			Object[] class1 = new Object[nbclasses1];
			Object[] class2 = new Object[nbclasses2];
			int i = 0;
			int j = 0;
			this.extractqqIfElseInnerForMethodAA(i, j, class1, class2);
			this.extractqqIfElseInnerForMethodAB(i, j, class1, class2, matrix, nbclasses1, nbclasses2);

			// Pass it to the algorithm
			int[][] result = callHungarianMethod( matrix, nbclasses1, nbclasses2 );
			// Extract the result
			this.extractqqIfElseInnerForMethodAC(i, result, class1, class2, threshold);
		}
	}

	/**
	 *
	 * @param i
	 * @param class1
	 * @param class2
	 */
	private void extractqqIfElseInnerForMethodAA(int i, int j, Object[] class1, Object[] class2){
		for ( Object ob : ontology1().getClasses() ) {
			class1[i++] = ob;
		}
		for ( Object ob : ontology2().getClasses() ) {
			class2[j++] = ob;
		}
	}

	/**
	 *
	 * @param i
	 * @param j
	 * @param class1
	 * @param class2
	 * @param matrix
	 * @param nbclasses1
	 * @param nbclasses2
	 */
	private void extractqqIfElseInnerForMethodAB(int i, int j, Object[] class1, Object[] class2, double[][] matrix, int nbclasses1, int nbclasses2){
		for( i = 0; i < nbclasses1; i++ ){
			for( j = 0; j < nbclasses2; j++ ){
				if ( sim.getSimilarity() ) matrix[i][j] = sim.getClassSimilarity(class1[i],class2[j]);
				else matrix[i][j] = 1. - sim.getClassSimilarity(class1[i],class2[j]);
			}
		}
	}

	/**
	 *
	 * @param i
	 * @param result
	 * @param class1
	 * @param class2
	 * @param threshold
	 */
	private void extractqqIfElseInnerForMethodAC(int i, int[][] result, Object[] class1, Object[] class2, double threshold){
		for( i=0; i < result.length ; i++ ){
			// The matrix has been destroyed
			double val;
			if ( sim.getSimilarity() ) val = sim.getClassSimilarity(class1[result[i][0]],class2[result[i][1]]);
			else val = 1 - sim.getClassSimilarity(class1[result[i][0]],class2[result[i][1]]);
			// JE: here using strict-> is a very good idea.
			// it means that correspondences with 0. similarity
			// will be excluded from the best match.
			if( val > threshold ){
				addCell( new ObjectCell( (String)null, class1[result[i][0]], class2[result[i][1]], BasicRelation.createRelation("="), val ) );
			}
		}
	}

	/**
	 *
	 * @param nbprop1
	 * @param nbprop2
	 * @param threshold
	 */
	private void extractqqIfElseMethodB(int nbprop1, int nbprop2, double threshold){
		if ( nbprop1 != 0 && nbprop2 != 0 ) {
			double[][] matrix = new double[nbprop1][nbprop2];
			Object[] prop1 = new Object[nbprop1];
			Object[] prop2 = new Object[nbprop2];
			int i = 0;
			ConcatenatedIterator pit1 = new
					ConcatenatedIterator(ontology1().getObjectProperties().iterator(),
					ontology1().getDataProperties().iterator());
			for ( Object ob: pit1 ) prop1[i++] = ob;
			int j = 0;
			ConcatenatedIterator pit2 = new
					ConcatenatedIterator(ontology2().getObjectProperties().iterator(),
					ontology2().getDataProperties().iterator());
			for ( Object ob: pit2 ) prop2[j++] = ob;
			this.extractqqIfElseInnerForMethodBA(i, j, matrix, nbprop1, nbprop2, prop1, prop2);

			// Pass it to the algorithm
			int[][] result = callHungarianMethod( matrix, nbprop1, nbprop2 );
			// Extract the result
			this.extractqqIfElseInnerForMethodBB( i, result, prop1, prop2, threshold);

		}
	}

	/**
	 *
	 * @param i
	 * @param j
	 * @param matrix
	 * @param nbprop1
	 * @param nbprop2
	 * @param prop1
	 * @param prop2
	 */
	private void extractqqIfElseInnerForMethodBA(int i, int j, int[][] matrix, int nbprop1, int nbprop2, Object[] prop1, Object[] prop2){
		for( i = 0; i < nbprop1; i++ ){
			for( j = 0; j < nbprop2; j++ ){
				if ( sim.getSimilarity() ) matrix[i][j] = sim.getPropertySimilarity(prop1[i],prop2[j]);
				else matrix[i][j] = 1. - sim.getPropertySimilarity(prop1[i],prop2[j]);
			}
		}
	}

	private void extractqqIfElseInnerForMethodBB(int i, int[][] result, Object[] prop1, Object[] prop2, double threshold){
		for( i=0; i < result.length ; i++ ){
			// The matrix has been destroyed
			double val;
			if ( sim.getSimilarity() ) val = sim.getPropertySimilarity(prop1[result[i][0]],prop2[result[i][1]]);
			else val = 1 - sim.getPropertySimilarity(prop1[result[i][0]],prop2[result[i][1]]);
			// JE: here using strict-> is a very good idea.
			// it means that alignments with 0. similarity
			// will be excluded from the best match.
			if( val > threshold ){
				addCell( new ObjectCell( (String)null, prop1[result[i][0]], prop2[result[i][1]], BasicRelation.createRelation("="), val ) );
			}
		}
	}

	/**
	 *
	 * @param params
	 * @param threshold
	 */
	private void extractqqIfElseMethodC(Properties params, double threshold){
		if (  params.getProperty("noinst") == null ){
			try {
				// Create individual lists
				Object[] ind1 = new Object[ontology1().nbIndividuals()];
				Object[] ind2 = new Object[ontology2().nbIndividuals()];
				int nbind1 = 0;
				int nbind2 = 0;
				this.extractqqIfElseInnerForMethodCA(ind1, nbind1);
				this.extractqqIfElseInnerForMethodCB(ind2, nbind2);
				this.extractqqIfElseInnerIfCC( nbind1, nbind2, ind1, ind2, threshold);

			} catch (AlignmentException alex) { alex.printStackTrace(); //}
			} catch (OntowrapException owex) { owex.printStackTrace(); }
		}
	}

	/**
	 *
	 * @param ind2
	 * @param nbind2
	 */
	private void extractqqIfElseInnerForMethodCA(Object[] ind2, int nbind2){
		for( Object ob : ontology2().getIndividuals() ){
			// We suppress anonymous individuals... this is not legitimate
			if ( ontology2().getEntityURI( ob ) != null ) {
				ind2[nbind2++] = ob;
			}
		}
	}

	/**
	 *
	 * @param ind1
	 * @param nbind1
	 */
	private void extractqqIfElseInnerForMethodCB(Object[] ind1, int nbind1){
		for( Object ob : ontology1().getIndividuals() ){
			// We suppress anonymous individuals... this is not legitimate
			if ( ontology1().getEntityURI( ob ) != null ) {
				ind1[nbind1++] = ob;
			}
		}
	}

	/**
	 *
	 * @param nbind1
	 * @param nbind2
	 * @param ind1
	 * @param ind2
	 * @param threshold
	 */
	private void extractqqIfElseInnerIfCC(int nbind1, int nbind2, Object[] ind1, Object[] ind2, double threshold){
		if ( nbind1 != 0 && nbind2 != 0 ) {
			double[][] matrix = new double[nbind1][nbind2];
			int i = 0;
			int j = 0;
			this.extractqqIfElseInnerForMethodCCA(matrix, nbind1, nbind2, i, j, ind1, ind2);

			// Pass it to the algorithm
			int[][] result = callHungarianMethod( matrix, nbind1, nbind2 );
			// Extract the result
			this.extractqqIfElseInnerForMethodCCB(int i, int[][] result, Object[] ind1, Object[] ind2, double threshold);

		}
	}

	/**
	 *
	 * @param matrix
	 * @param nbind1
	 * @param nbind2
	 * @param i
	 * @param j
	 * @param ind1
	 * @param ind2
	 */
	private void extractqqIfElseInnerForMethodCCA(double[][] matrix, int nbind1, int nbind2, int i, int j, Object[] ind1, Object[] ind2){
		for( i=0; i < nbind1; i++ ){
			for( j=0; j < nbind2; j++ ){
				if ( sim.getSimilarity() ) matrix[i][j] = sim.getIndividualSimilarity(ind1[i],ind2[j]);
				else matrix[i][j] = 1 - sim.getIndividualSimilarity(ind1[i],ind2[j]);
			}
		}
	}

	/**
	 *
	 * @param i
	 * @param result
	 * @param ind1
	 * @param ind2
	 * @param threshold
	 */
	private void extractqqIfElseInnerForMethodCCB(int i, int[][] result, Object[] ind1, Object[] ind2, double threshold){
		for( i=0; i < result.length ; i++ ){
			// The matrix has been destroyed
			double val;
			if ( sim.getSimilarity() ) val = sim.getIndividualSimilarity(ind1[result[i][0]],ind2[result[i][1]]);
			else val = 1 - sim.getIndividualSimilarity(ind1[result[i][0]],ind2[result[i][1]]);
			// JE: here using strict-> is a very good idea.
			// it means that alignments with 0. similarity
			// will be excluded from the best match.
			if( val > threshold ){
				addCell( new ObjectCell( (String)null, ind1[result[i][0]], ind2[result[i][1]], BasicRelation.createRelation("="), val ) );
			}
		}
	}


    /**
     * Greedy algorithm:
     * 1) dump the part of the matrix distance above threshold in a sorted set
     * 2) traverse the sorted set and each time a correspondence involving two
     *    entities that have no correspondence is encountered, add it to the 
     *    alignment.
     * Complexity: O(n^2.logn)
     * Pitfall: no global optimality is warranted, nor stable marriage
     * for instance if there is the following matrix:
     * (a,a')=1., (a,b')=.9, (b,a')=.9, (b,b')=.1
     * This algorithm will select the first and last correspondances of
     * overall similarity 1.1, while the optimum is the second solution
     * with overall of 1.8.
     */
    @SuppressWarnings("unchecked") //ConcatenatedIterator
    public Alignment extractqqgreedy( double threshold, Properties params) {
	double val = 0;
	//TreeSet could be replaced by something else
	//The comparator must always tell that things are different!
	SortedSet<Cell> cellSet = new TreeSet<Cell>(
			    new Comparator<Cell>() {
				public int compare( Cell o1, Cell o2 )
				    throws ClassCastException{
				    try {
					//System.err.println(o1.getObject1()+" -- "+o1.getObject2()+" // "+o2.getObject1()+" -- "+o2.getObject2());
					if ( o1.getStrength() > o2.getStrength() ){
					    return -1;
					} else if ( o1.getStrength() < o2.getStrength() ){
					    return 1;
					} else if ( ontology1().getEntityName( o1.getObject1() ) == null
						    || ontology2().getEntityName( o2.getObject1() ) == null ) {
					    return -1;
					} else if ( ontology1().getEntityName( o1.getObject1()).compareTo( ontology2().getEntityName( o2.getObject1() ) ) > 0 ) {
					    return -1;
					} else if ( ontology1().getEntityName( o1.getObject1()).compareTo( ontology2().getEntityName( o2.getObject1() ) ) < 0 ) {
					    return 1;
					} else if ( ontology1().getEntityName( o1.getObject2() ) == null
						    || ontology2().getEntityName( o2.getObject2() ) == null ) {
					    return -1;
					} else if ( ontology1().getEntityName( o1.getObject2()).compareTo( ontology2().getEntityName( o2.getObject2() ) ) > 0 ) {
					    return -1;
					// Assume they have different names
					} else { return 1; }
				    } catch ( OntowrapException e) { 
					e.printStackTrace(); return 0;}
				}
			    }
			    );
      try {
	  // Get all the matrix above threshold in the SortedSet
	  // Plus a map from the objects to the cells
	  // O(n^2.log n)
	  // for classes
	  for ( Object ent1: ontology1().getClasses() ) {
	      for ( Object ent2: ontology2().getClasses() ) {
		  if ( sim.getSimilarity() ) val = sim.getClassSimilarity( ent1, ent2 );
		  else val = 1 - sim.getClassSimilarity( ent1, ent2 );
		  if ( val > threshold ){
		      cellSet.add( new ObjectCell( (String)null, ent1, ent2, BasicRelation.createRelation("="), val ) );
		  }
	      }
	  }
	  // for properties
	  ConcatenatedIterator pit1 = new 
	      ConcatenatedIterator(ontology1().getObjectProperties().iterator(),
				   ontology1().getDataProperties().iterator());
	  for ( Object ent1: pit1 ) {
	      ConcatenatedIterator pit2 = new 
		  ConcatenatedIterator(ontology2().getObjectProperties().iterator(),
					ontology2().getDataProperties().iterator());
	      for ( Object ent2: pit2 ) {
		  if ( sim.getSimilarity() ) val = sim.getPropertySimilarity( ent1, ent2 );
		  else val = 1 - sim.getPropertySimilarity( ent1, ent2 );
		  if ( val > threshold ){
		      cellSet.add( new ObjectCell( (String)null, ent1, ent2, BasicRelation.createRelation("="), val ) );
		  }
	      }
	  }
	  // for individuals
	  if (  params.getProperty("noinst") == null ){
	      for( Object ent1: ontology1().getIndividuals() ) {
		  if ( ontology1().getEntityURI( ent1 ) != null ) {

		      for( Object ent2: ontology2().getIndividuals() ) {
			  if ( ontology2().getEntityURI( ent2 ) != null ) {
			      if ( sim.getSimilarity() ) val = sim.getIndividualSimilarity( ent1, ent2 );
			      else val = 1 - sim.getIndividualSimilarity( ent1, ent2 );
			      if ( val > threshold ){
				  cellSet.add( new ObjectCell( (String)null, ent1, ent2, BasicRelation.createRelation("="), val ) );
			      }
			  }
		      }
		  }
	      }
	  }

	  // O(n^2)
	  for( Cell cell : cellSet ){
	      Object ent1 = cell.getObject1();
	      Object ent2 = cell.getObject2();
	      if ( (getAlignCells1( ent1 ) == null) && (getAlignCells2( ent2 ) == null) ){
		  // The cell is directly added!
		  addCell( cell );
	      }
	  };

      } catch (AlignmentException alex) {
	  alex.printStackTrace();
      } catch (OntowrapException owex) {
	  owex.printStackTrace();
      };
      return((Alignment)this);
    }

}
