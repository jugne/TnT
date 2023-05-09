package tnt.util;

import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Taxon;

public class MemoryFriendlyAlignment extends FilteredAlignment {

    @Override
    public void finalize(){
        sequences.clear();
        patternWeight = null;
        patternIndex = null;
        sitePatterns = null;
        siteWeights = null;
        taxaNames.clear();
        m_dataType = null;
        stateCounts.clear();
        sequenceInput.get().clear();
        sequenceInput.defaultValue.clear();
        taxonSetInput.setValue(null, new Taxon());
        taxonSetInput.defaultValue.taxonsetInput.get().clear();
    }
}
