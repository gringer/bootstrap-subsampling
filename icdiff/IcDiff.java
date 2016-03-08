//package icdiff;

import java.util.zip.*;
import java.util.Vector;
import java.util.HashMap;
import java.io.*;

/**
 * @author David Hall (gringer)
 * @version 0.1
 *
 * Created 2009-Jun-04
 *
 * This program calculates information content difference for a
 * chromosome. The default data format is assumed to be one marker per
 * line, containing (possibly) a marker location indicator following
 * the marker.
*/
public class IcDiff {
    /**
     * @param args
     *            List of files to process
     */
    static int POS_GT = 2;

    static int POS_LOCATION = 1;

    static int POS_MARKER = 0;

    static boolean DEBUG = false;

    public static double icDiff(Vector<StringBuffer> haplotypes, int location) {
        return IcDiff.icDiff(haplotypes, location, true);
    }

    public static double icDiff(Vector<StringBuffer> haplotypes, int location, boolean constrain){
        if(haplotypes.size() == 0){
            return 0.0;
        }
        double icTarget = Math.log((double)haplotypes.size()) / Math.log(2.0);
        double icDiffSum = 0.0;
        int numMarkers = haplotypes.get(0).length();
        for(int leftValue = 0; leftValue < 2; leftValue ++){
            boolean goLeft = (leftValue == 0); // having a boolean in the for didn't work
            if (IcDiff.DEBUG) {
                System.out.println("Going left is " + goLeft);
            }
            double icValue = 0.0;
            double icMax = 0.0;
            int locLeft = location;
            int locRight = location;
            int numDifferent = 0;
            while(numDifferent < haplotypes.size()){
                if(constrain){
                    // constrain information content to maximum possible with
                    // given number of haplotypes
                    icMax = Math.min(((locRight - locLeft) + 1), icTarget);
                } else {
                    icMax = (locRight - locLeft) + 1;
                }
                HashMap<String,Integer> probHap = new HashMap<String,Integer>();
                // This loop is the following converted to Java
                /* Vector<Double> probHap = table(apply(haplotypes[,locLeft:locRight],
                1, paste, collapse = "")) / numHaplotypes; */
                for(int i = 0; i < haplotypes.size(); i++){
                    String key = haplotypes.get(i).substring(locLeft, locRight+1);
                    if(probHap.containsKey(key)){
                        probHap.put(key,probHap.get(key) + 1);
                    } else {
                        probHap.put(key,1);
                    }
                }
                // This loop is the following converted to Java
                /* icValue = sum(-probHap * Math.log(probHap)/Math.log(2)); */
                icValue = 0.0;
                for(String key : probHap.keySet()){
                    double probValue = (double)probHap.get(key) / haplotypes.size();
                    if (IcDiff.DEBUG) {
                        System.out.print(key + " occurrs " + probHap.get(key) + " times");
                        System.out.println("(" + probValue + ")");
                    }
                    probValue = -probValue * (Math.log(probValue) / Math.log(2));
                    icValue += probValue;
                }
                numDifferent = probHap.size();
                if (IcDiff.DEBUG) {
                    System.out.println("Value is now " + icValue + " (target = " + icTarget + ")");
                }
                icDiffSum += ((icMax - icValue) / icMax);
                // Extra logic to make sure it doesn't duplicate results by
                // hanging in the same place for two iterations. These two
                // 'goLeft' cases could be combined, but it makes the code less
                // clear.
                if (IcDiff.DEBUG) {
                    System.out.println("Go left state is " + goLeft);
                }
                if(goLeft){
                    goLeft = !goLeft; // i.e. goLeft = FALSE
                    if(locLeft > 0){
                        locLeft = Math.max(0, locLeft - 1);
                    } else {
                        locRight = Math.min(numMarkers - 1, locRight + 1);
                    }
                } else {
                    goLeft = !goLeft; // i.e. goLeft = TRUE
                    if(locRight < (numMarkers - 1)){
                        locRight = Math.min(numMarkers - 1, locRight + 1);
                    } else {
                        locLeft = Math.max(0, locLeft - 1);
                    }
                }
                if((locLeft == 1) && (locRight == numMarkers)){
                    // if the entire chromosome is covered, do no more
                    icValue = icTarget;
                }
            }
        }
        return(icDiffSum / 2);
    }

    public static int readToken(StreamTokenizer st){
        int tokenID = StreamTokenizer.TT_EOF;
        try {
            tokenID = st.nextToken();
        } catch (IOException e) {
            /*
             * workaround a Java bug caused by GZIP files being
             * >2GB
             * http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=4262583
             * ... Of course, there's a problem here if the GZIP
             * file actually <em>does</em> have a corrupt
             * trailer.
             */
            if (e.getMessage().contains("Corrupt GZIP trailer")) {
                tokenID = StreamTokenizer.TT_EOF;
            } else {
                System.err.println("Error in reading token from file:");
                System.err.println(e.getMessage());
                e.printStackTrace();
                tokenID = StreamTokenizer.TT_EOF;
            }
        }
        return tokenID;
    }

    public static void flushLine(StreamTokenizer st){
        // clears a StreamTokenizer to the end of the line
        int tokenID = -4;
        while((tokenID != StreamTokenizer.TT_EOL) && (tokenID != StreamTokenizer.TT_EOF)){
            // pass over remainder of line
            tokenID = IcDiff.readToken(st);
        }
    }

    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        int linesRead = 0;
        // boolean firstTime = true;
        boolean fastphase = false; // intially assume not in fastphase format
        boolean fastphaseBegin = false;
        if (IcDiff.DEBUG) {
            System.out.println("Reading in files");
        }
        StreamTokenizer inputScanner = null;
        /* Assume 100 haplotypes initially. This will not change after the first line */
        Vector<StringBuffer> haplotypes = new Vector<StringBuffer>(100);
        /* Assume 30000 markers initially */
        Vector<String> markers = new Vector<String>(30000);
        for (int argPos = 0; argPos < args.length; argPos++) {
            haplotypes.setSize(0); // reset size of haplotypes
            int tokenCounter = 0;
            linesRead = 0;
            System.err.println("Reading " + args[argPos]);
            try {
                inputScanner = new StreamTokenizer(new InputStreamReader(
                        new GZIPInputStream(new FileInputStream(args[argPos]))));
            } catch (FileNotFoundException e) {
                System.out.println("Error: file not found '" + args[argPos]
                        + "'");
                System.exit(1);
            } catch (IOException e1) {
                if (e1.getMessage().contains("Not in GZIP format")) {
                    try {
                        inputScanner = new StreamTokenizer(
                                new InputStreamReader(new FileInputStream(
                                        args[argPos])));
                    } catch (FileNotFoundException e) {
                        System.out.println("Error: file not found '"
                                + args[argPos] + "'");
                        System.exit(1);
                    }
                } else {
                    e1.printStackTrace();
                }
            }
            /* Consider end-of-line to be a special token */
            inputScanner.eolIsSignificant(true);
            /*
             * Consider non-numeric printable ASCII characters to be word/number
             * characters
             */
            inputScanner.wordChars(33, 126);
            if (inputScanner != null) {
                int tokenPos = 0;
                // String curMarker = null;
                int tTokenID = StreamTokenizer.TT_EOL;
                while (tTokenID != StreamTokenizer.TT_EOF) {
                    /*
                     * The general process is to read the input files, line
                     * by line, and concatenate genotypes into single
                     * haplotypes for subsequent calculation. Data are
                     * stored in a vector of Strings.
                     */
                    tTokenID = IcDiff.readToken(inputScanner);
                    /*
                     * Number: probably marker location or similar
                     */
                    if ((!fastphase) && (tTokenID == StreamTokenizer.TT_NUMBER)) {
                        if (IcDiff.DEBUG) {
                            System.out.println("Read number [" + tokenPos
                                               + "]: " + inputScanner.nval);
                        }
                        tokenCounter++;
                        tokenPos++;
                    }
                    /*
                     * Word: Some combination of alphanumeric characters,
                     * hyphens, and other odd symbols. This is treated as
                     * either a marker name, or a genotype, depending on the
                     * field/column number.
                     */
                    if (tTokenID == StreamTokenizer.TT_WORD) {
                        if (IcDiff.DEBUG) {
                            System.err.println("Read word [" + tokenPos
                                               + "]: " + inputScanner.sval);
                        }
                        if (tokenPos == IcDiff.POS_MARKER) {
                            if(fastphase){
                                if(fastphaseBegin && inputScanner.sval.equals("#")){
                                    // found an ID line
                                    IcDiff.flushLine(inputScanner);
                                    tokenPos = 0;
                                    tTokenID = -4;
                                    StringBuffer tmpHap = new StringBuffer("");
                                    while((tTokenID != StreamTokenizer.TT_EOL) && (tTokenID != StreamTokenizer.TT_EOF)){
                                        // read in alleles for first haplotype
                                        tTokenID = IcDiff.readToken(inputScanner);
                                        if(tTokenID != StreamTokenizer.TT_EOL){
                                            if(markers.size() <= tokenPos){
                                                markers.add(Integer.toString(tokenPos+1)); // R code is 1-based
                                            }
                                            if(tTokenID == StreamTokenizer.TT_NUMBER){
                                                tmpHap.append((int)inputScanner.nval);
                                            } else {
                                                tmpHap.append(inputScanner.sval);
                                            }
                                            tokenPos++;
                                        }
                                    }
                                    haplotypes.add(tmpHap);
                                    tokenPos = 0;
                                    tTokenID = -4;
                                    tmpHap = new StringBuffer("");
                                    while((tTokenID != StreamTokenizer.TT_EOL) && (tTokenID != StreamTokenizer.TT_EOF)){
                                        // read in line for second haplotype
                                        tTokenID = IcDiff.readToken(inputScanner);
                                        if(tTokenID != StreamTokenizer.TT_EOL){
                                            if(tTokenID == StreamTokenizer.TT_NUMBER){
                                                tmpHap.append((int)inputScanner.nval);
                                            } else {
                                                tmpHap.append(inputScanner.sval);
                                            }
                                            tokenPos++;
                                        }
                                    }
                                    haplotypes.add(tmpHap);
                                    tokenPos = 0;
                                } else if(inputScanner.sval.equals("BEGIN")){
                                    tTokenID = IcDiff.readToken(inputScanner);
                                    if(inputScanner.sval.equals("GENOTYPES")){
                                        fastphaseBegin = true;
                                    }
                                    IcDiff.flushLine(inputScanner);
                                } else if(inputScanner.sval.equals("END")){
                                    fastphaseBegin = false;
                                    IcDiff.flushLine(inputScanner);
                                } else {
                                    IcDiff.flushLine(inputScanner);
                                }
                            } else {
                                if(inputScanner.sval.substring(0,2).equals("**")){
                                    // assume fastphase output
                                    System.err.println("Assuming fastphase output...");
                                    fastphase = true;
                                } else {
                                    markers.add(inputScanner.sval);
                                }
                            }
                        } else if ((!fastphase) && (tokenPos >= IcDiff.POS_GT)) {
                            int indPos = tokenPos - IcDiff.POS_GT;
                            if (indPos >= haplotypes.size()) {
                                if (IcDiff.DEBUG) {
                                    System.err
                                        .println("New person needed at position "
                                                 + indPos);
                                }
                                haplotypes.add(new StringBuffer(
                                                                inputScanner.sval));
                            } else {
                                StringBuffer tmpHap = haplotypes
                                    .get(indPos);
                                tmpHap.append(inputScanner.sval);
                            }
                        }
                        tokenCounter++;
                        if(!fastphase){
                            // don't update position counter for fastphase files
                            // -- they use a different input parsing method
                            tokenPos++;
                        }
                    }
                    if (tTokenID == StreamTokenizer.TT_EOL) {
                        if (IcDiff.DEBUG) {
                            System.err.println("[New Line]");
                        }
                        tokenPos = 0;
                        linesRead++;
                    }
                }
                // remove similar haplotypes
                Vector<Integer> removeHaps = new Vector<Integer>(haplotypes.size());
                // Minimum number of differences for haplotypes to be considered different
                // (This allows for some genotyping error)
                int minDifferences = 20;
                System.err.print("Checking for similar haplotypes (< " + minDifferences +
                                   " differences)... ");
                for (int hapPos1 = 0; hapPos1 < haplotypes.size(); hapPos1++) {
                    String compare1 = haplotypes.get(hapPos1).toString();
                    int countDifferences = 0;
                    if(removeHaps.contains(hapPos1)){
                        continue; // already removed, so don't compare
                    }
                    for (int hapPos2 = hapPos1+1; hapPos2 < haplotypes.size(); hapPos2++) {
                        String compare2 = haplotypes.get(hapPos2).toString();
                        for(int i = 0; i < Math.min(compare1.length(), compare2.length()); i++){
                            if(compare1.charAt(i) != compare2.charAt(i)){
                                countDifferences++;
                                if(countDifferences >= minDifferences){
                                    break; // no need to look for more differences
                                }
                            }
                        }
                        if(countDifferences < minDifferences){
                            removeHaps.add(hapPos2); // add to list of haplotypes to remove
                        }
                    }
                }
                System.err.println("done");
                if(removeHaps.size() > 0){
                    System.err.print("About to remove " + removeHaps.size() + " haplotypes... ");
                }
                for(int i = removeHaps.size() - 1; i >= 0; i--){
                    // done in reverse to preserve location information
                    haplotypes.remove(removeHaps.get(i).intValue());
                }
                if(removeHaps.size() > 0){
                    System.err.println("done, " + haplotypes.size() + " haplotypes remain");
                }
                if (IcDiff.DEBUG) {
                    for (int hapPos = 0; hapPos < haplotypes.size(); hapPos++) {
                        System.out.println(haplotypes.get(hapPos));
                    }
                }
                // now do icDiff calculations
                for (int markerPos = 0; markerPos < markers.size(); markerPos++) {
                    //messy hack for 7 significant figures
                    double markerDiff = IcDiff.icDiff(haplotypes, markerPos);
                    int magnitude = (int)Math.log10(markerDiff);
                    int dp = 7 - magnitude - 1;
                    System.out.printf("%s %."+dp+"f", markers.get(markerPos),
                                      markerDiff);
                    System.out.println();
                }
                System.err
                    .println("Completed in "
                             + ((System.currentTimeMillis() - startTime) / 1000.0)
                             + " seconds");
            }
        }
    }
}
