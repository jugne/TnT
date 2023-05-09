package tnt.util;



import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MemoryFriendlyAlignment extends FilteredAlignment {

public void freeMemory(){
    sequences.clear();
    patternWeight = null;
    patternIndex = null;
    sitePatterns = null;
    siteWeights = null;
    taxaNames.clear();
    m_dataType = null;
//    counts.clear();
    stateCounts.clear();
    sequenceInput.get().clear();
    sequenceInput.defaultValue.clear();
    taxonSetInput.setValue(null, new Taxon());
//    taxonSetInput.set(null);
    taxonSetInput.defaultValue.taxonsetInput.get().clear();
}




}
