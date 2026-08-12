@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')
@Grab(group='org.codehaus.gpars', module='gpars', version='1.1.0')

import org.openrdf.model.URI
import org.openrdf.model.vocabulary.*
import slib.sglib.io.loader.*
import slib.sglib.algo.graph.utils.*
import slib.graph.algo.extraction.utils.*
import slib.graph.algo.extraction.rvf.instances.*
import slib.graph.algo.extraction.rvf.instances.impl.*
import slib.graph.model.graph.*
import slib.graph.model.repo.*
import slib.graph.model.impl.graph.memory.*
import slib.graph.model.impl.graph.elements.*
import slib.graph.model.impl.repo.*
import slib.graph.io.conf.*
import slib.graph.io.util.*
import slib.graph.io.loader.*
import groovyx.gpars.GParsPool

System.setProperty("jdk.xml.entityExpsansionLimit", "0")
System.setProperty("jdk.xml.totalEntitySizeLimit", "0")

def factory = URIFactoryMemory.getSingleton()
def inputPath  = args[0]  // TSV: Protein<tab>go_id<tab>go_ic (with header)
def resSimPath = args[1]  // output similarity matrix
def owlPath    = args[2]  // path to go-basic.owl (or equivalent OWL file)

def goIdToUri = { String goId ->
    def parts = goId.split(':')
    factory.getURI("http://purl.obolibrary.org/obo/${parts[0]}_${parts[1]}")
}

// Parse input: build protein->annotations and go_uri->ic maps
def proteinAnnotations = new LinkedHashMap()   // preserves insertion order
def icMap = [:]

new File(inputPath).eachLine { line ->
    if (line.startsWith("Protein")) return
    def parts = line.trim().split('\t')
    if (parts.size() < 3) return
    def protein = parts[0]
    def goUri   = goIdToUri(parts[1])
    def ic      = parts[2].toDouble()
    if (!proteinAnnotations.containsKey(protein)) proteinAnnotations[protein] = new LinkedHashSet()
    proteinAnnotations[protein].add(goUri)
    icMap[goUri] = ic
}

def proteins = proteinAnnotations.keySet().toList()
println "Proteins: ${proteins}"

// Load GO ontology (needed for ancestor traversal to find MICA)
URI graphUri = factory.getURI("http://purl.obolibrary.org/obo/")
def graph = new GraphMemory(graphUri)
GDataConf goConf = new GDataConf(GFormat.RDF_XML, owlPath)
GraphLoaderGeneric.populate(goConf, graph)

URI virtualRoot = factory.getURI("http://purl.obolibrary.org/obo/virtualRoot")
graph.addV(virtualRoot)
GAction rooting = new GAction(GActionType.REROOTING)
rooting.addParameter("root_uri", virtualRoot.stringValue())
GraphActionExecutor.applyAction(factory, rooting, graph)

// Cache inclusive ancestors (term + all superclasses) via RDFS.SUBCLASSOF traversal
def ancestorCache = [:]
def getAncestorsInc = { URI uri ->
    if (ancestorCache.containsKey(uri)) return ancestorCache[uri]
    def visited = new LinkedHashSet()
    def queue   = new LinkedList([uri])
    while (!queue.isEmpty()) {
        URI v = queue.poll()
        if (visited.contains(v)) continue
        visited.add(v)
        graph.getE(v).each { edge ->
            if (edge.getSource() == v && edge.getPredicate() == RDFS.SUBCLASSOF)
                queue.add(edge.getTarget())
        }
    }
    ancestorCache[uri] = visited
    return visited
}

// Resnik pairwise: IC of the Most Informative Common Ancestor (MICA)
def resnikSim = { URI a, URI b ->
    def common = getAncestorsInc(a).intersect(getAncestorsInc(b))
    if (common.isEmpty()) return 0.0
    def mica_ic = common.collect { icMap.getOrDefault(it, 0.0) }.max() ?: 0.0
    def ic_a = icMap.getOrDefault(a, 0.0)
    def ic_b = icMap.getOrDefault(b, 0.0)
    return (ic_a + ic_b) > 0 ? 2.0 * mica_ic / (ic_a + ic_b) : 0.0
}

// BMA (Best Match Average) groupwise similarity
def bmaSim = { Set annotsA, Set annotsB ->
    if (!annotsA || !annotsB) return 0.0
    double sum = 0.0
    annotsA.each { a -> sum += annotsB.collect { b -> resnikSim(a, b) }.max() }
    annotsB.each { b -> sum += annotsA.collect { a -> resnikSim(a, b) }.max() }
    return sum / (annotsA.size() + annotsB.size())
}

// Precompute ancestor sets for all GO terms in our data (before parallel section)
icMap.keySet().each { uri -> getAncestorsInc(uri) }

// Compute N×N similarity matrix in parallel
def n = proteins.size()
def result = new double[n][n]

GParsPool.withPool {
    (0..<n).eachParallel { i ->
        for (int j = 0; j < n; j++) {
            result[i][j] = bmaSim(proteinAnnotations[proteins[i]], proteinAnnotations[proteins[j]])
        }
    }
}

// Write labelled matrix
def out = new PrintWriter(new BufferedWriter(new FileWriter(resSimPath)))
out.println("\t" + proteins.join("\t"))
for (int i = 0; i < n; i++) {
    def row = [proteins[i]] + (0..<n).collect { j -> String.format("%.4f", result[i][j]) }
    out.println(row.join("\t"))
}
out.flush()
out.close()

println "Done. Matrix written to ${resSimPath}"
