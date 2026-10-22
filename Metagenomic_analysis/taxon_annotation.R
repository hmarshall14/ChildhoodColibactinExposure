###############################################################################
# taxon_annotation.R
#
# Shared annotation module for the crypt-WGS metagenome scripts.
# source() this from 01_coverage_diagnostics.R and 02_association_testing.R.
#
# Produces THREE INDEPENDENT annotation columns per taxon, deliberately not
# collapsed into one verdict:
#
#   blacklist_status  contaminant / conflicted / clean
#                     "is this a documented reagent, kit, skin or lab organism?"
#   allowlist_status  gut / oral / unknown / unnamed_sgb
#                     "is there positive evidence this colonises the human gut?"
#   evidence_score    data-driven, from your own counts (depth correlation,
#                     blanks, coverage breadth, prevalence)
#
# Read them together. Disagreement is informative: a taxon that is on the
# blacklist AND in the gut catalogue (E. coli, K. pneumoniae) cannot be
# resolved from a taxon table at all and must be settled at read level.
#
# SOURCES for the curated lists (see each block):
#   Salter et al. 2014, BMC Biology 12:87           (kit/reagent contaminants)
#   Eisenhofer et al. 2019, Trends Microbiol 27:105 (consolidated contaminant list)
#   Davis et al. 2018, Microbiome 6:226             (decontam)
#   Almeida et al. 2021, Nat Biotechnol 39:105      (UHGG gut genome catalogue)
#   Gihawi et al. 2023, mBio / Sci Rep              (cautionary: tissue-WGS artefacts)
#
# The lists below are a curated working set, NOT exhaustive. For a proper
# allowlist, download UHGG and pass it to annotate_taxa(gut_catalogue = ...).
###############################################################################

## ===========================================================================
## A. GENUS-LEVEL CONTAMINANTS
##    Genera in which essentially every species seen in low-biomass work is
##    a reagent/environmental organism. Safe to flag wholesale.
## ===========================================================================

GENUS_CONTAM <- c(
  # rhizosphere / plant-associated -- extremely common in kit and water blanks
  "Bradyrhizobium","Rhizobium","Neorhizobium","Pararhizobium","Allorhizobium",
  "Agrobacterium","Mesorhizobium","Sinorhizobium","Ensifer","Phyllobacterium",
  "Devosia","Aureimonas","Methylobacterium","Methylorubrum","Bosea",
  "Nitrobacter","Rhodopseudomonas","Rhodoplanes","Afipia","Ochrobactrum",
  "Brucella","Xanthomonas","Pseudoxanthomonas","Lysobacter","Dyella",
  # beta/gamma-proteobacterial water organisms
  "Ralstonia","Cupriavidus","Delftia","Variovorax","Acidovorax","Comamonas",
  "Hydrogenophaga","Polaromonas","Curvibacter","Diaphorobacter","Limnobacter",
  "Undibacterium","Herbaspirillum","Janthinobacterium","Massilia","Duganella",
  "Pelomonas","Schlegelella","Tepidimonas","Ideonella","Rubrivivax",
  "Aquabacterium","Leptothrix","Methylophilus","Methylotenera","Perlucidibaca",
  # sphingomonads / caulobacters
  "Sphingomonas","Sphingobium","Sphingopyxis","Novosphingobium","Blastomonas",
  "Porphyrobacter","Erythrobacter","Caulobacter","Brevundimonas",
  "Asticcacaulis","Phenylobacterium","Phreatobacter","Craurococcus",
  "Roseomonas","Acidiphilium","Paracoccus","Rubellimicrobium","Chelatococcus",
  # actinobacterial skin / environment
  "Cutibacterium","Propionibacterium","Microbacterium","Leifsonia","Curtobacterium",
  "Kocuria","Dermacoccus","Brevibacterium","Rhodococcus","Nocardioides",
  "Williamsia","Gordonia","Mycolicibacterium","Nakamurella","Modestobacter",
  # misc environmental
  "Deinococcus","Thermus","Meiothermus","Pedobacter","Mucilaginibacter",
  "Hymenobacter","Spirosoma","Dyadobacter","Chitinophaga","Sediminibacterium",
  "Polynucleobacter","Halomonas","Marinobacter","Shewanella","Vibrio",
  "Alcanivorax","Cellvibrio","Luteimonas","Steroidobacter",
  # --- added after inspecting a full-NCBI Bracken profile of crypt WGS ------
  # Pseudomonas was split in 2021; the new genera inherit the contaminant
  # status of the species that moved into them (P. stutzeri -> Stutzerimonas).
  "Stutzerimonas","Halopseudomonas","Metapseudomonas","Atlantibacter",
  # halophilic archaea -- soda lakes and salterns. Cannot be colonic residents.
  # Large read counts here are a misassignment signal, not a biological one.
  "Natronomonas","Halopenitus","Haloferax","Halorubrum","Haloarcula",
  "Halobacterium","Natrialba","Natronococcus","Halomicrobium","Haloterrigena",
  "Natrinema","Halorhabdus","Haloquadratum","Halogeometricum","Halovivax",
  "Natronolimnobius","Salinibacter","Halalkalicoccus","Haloplanus",
  # soil Acidobacteria
  "Edaphobacter","Granulicella","Acidobacterium","Terriglobus","Bryocella",
  # soil/environmental Micrococcaceae and Microbacteriaceae
  "Arthrobacter","Pseudarthrobacter","Paenarthrobacter","Glutamicibacter",
  "Sanguibacter","Agromyces","Cellulomonas","Frigoribacterium","Plantibacter",
  "Clavibacter","Okibacterium","Rathayibacter","Kribbella","Streptacidiphilus",
  # environmental Flavobacteriaceae (hospital water, wet surfaces)
  "Elizabethkingia","Epilithonimonas","Empedobacter","Myroides","Wautersiella",
  "Sphingobacterium","Olivibacter",
  # soil/reagent spore formers
  "Paenibacillus","Brevibacillus","Lysinibacillus","Solibacillus","Sporosarcina",
  "Aneurinibacillus","Ureibacillus","Cohnella","Fictibacillus"
)

## Whole clades that cannot be human colonic residents. Matched against the
## FULL LINEAGE string (not the species name), so they catch the entire
## environmental long tail of a full-NCBI Kraken/Bracken database without
## needing every genus enumerated.
## NB: deliberately NOT blanket-excluding Euryarchaeota (Methanobrevibacter
## smithii is a genuine gut methanogen) or Verrucomicrobiota (Akkermansia).
ENV_CLADE_PATTERNS <- c(
  "c__Halobacteria", "o__Halobacteriales", "o__Haloferacales",
  "o__Natrialbales", "p__Halobacteriota", "p__Halobacterota",
  "p__Crenarchaeota", "p__Thaumarchaeota", "p__Nitrososphaerota",
  "p__Acidobacteri", "p__Deinococc", "p__Chloroflex", "p__Cyanobacteri",
  "p__Planctomycet", "p__Thermotog", "p__Aquific", "p__Nitrospir",
  "p__Gemmatimonadet", "p__Chlorobi", "p__Armatimonadet", "p__Caldiserica",
  "p__Dictyoglomi", "p__Elusimicrobi", "p__Rhodothermot", "p__Balneolaeot",
  "f__Halomonadaceae", "f__Rhodobacteraceae", "f__Sphingomonadaceae",
  "f__Xanthobacteraceae", "f__Nitrobacteraceae", "f__Beijerinckiaceae",
  "f__Acetobacteraceae", "f__Caulobacteraceae", "f__Hyphomicrobiaceae",
  "o__Chroococcales", "o__Nostocales", "o__Oscillatoriales",
  "g__Prochlorococcus", "g__Synechococcus"
)

## ===========================================================================
## B. SPECIES-LEVEL CONTAMINANTS
##    Named species repeatedly recovered from extraction kits, PCR reagents,
##    ultrapure water, flow cells and skin. Use when RANK = "s".
## ===========================================================================

SPECIES_CONTAM <- c(
  # --- skin / operator -----------------------------------------------------
  "Cutibacterium acnes","Cutibacterium granulosum","Cutibacterium avidum",
  "Propionibacterium acnes",
  "Staphylococcus epidermidis","Staphylococcus hominis","Staphylococcus capitis",
  "Staphylococcus warneri","Staphylococcus haemolyticus","Staphylococcus caprae",
  "Staphylococcus saccharolyticus","Staphylococcus pettenkoferi",
  "Corynebacterium tuberculostearicum","Corynebacterium kroppenstedtii",
  "Corynebacterium simulans","Corynebacterium afermentans",
  "Micrococcus luteus","Micrococcus lylae","Kocuria rhizophila",
  "Malassezia restricta","Malassezia globosa",
  "Moraxella osloensis","Dermabacter hominis",
  # --- water / reagent / kit ----------------------------------------------
  "Ralstonia pickettii","Ralstonia insidiosa","Ralstonia mannitolilytica",
  "Ralstonia solanacearum",
  "Delftia acidovorans","Delftia tsuruhatensis",
  "Stenotrophomonas maltophilia","Stenotrophomonas rhizophila",
  "Sphingomonas paucimobilis","Sphingomonas echinoides","Sphingomonas koreensis",
  "Sphingomonas melonis","Sphingomonas aerolata",
  "Methylobacterium radiotolerans","Methylobacterium extorquens",
  "Methylorubrum populi",
  "Herbaspirillum huttiense","Herbaspirillum seropedicae",
  "Afipia broomeae","Afipia felis",
  "Variovorax paradoxus","Acidovorax temperans","Acidovorax delafieldii",
  "Comamonas testosteroni","Comamonas kerstersii",
  "Burkholderia cepacia","Burkholderia fungorum","Burkholderia contaminans",
  "Paraburkholderia fungorum",
  "Caulobacter vibrioides","Caulobacter crescentus","Caulobacter segnis",
  "Pelomonas puraquae","Pelomonas saccharophila",
  "Brevundimonas diminuta","Brevundimonas vesicularis","Brevundimonas naejangsanensis",
  "Janthinobacterium lividum",
  "Achromobacter xylosoxidans","Achromobacter insuavis",
  "Rhizobium radiobacter","Agrobacterium tumefaciens","Rhizobium pusense",
  "Massilia timonae","Chryseobacterium indologenes","Chryseobacterium gleum",
  "Flavobacterium johnsoniae",
  "Cupriavidus metallidurans","Cupriavidus necator","Cupriavidus pauculus",
  "Paracoccus yeei","Paracoccus denitrificans",
  "Phyllobacterium myrsinacearum",
  "Deinococcus radiodurans",
  "Microbacterium oxydans","Microbacterium testaceum","Leifsonia xyli",
  "Bacillus subtilis","Bacillus licheniformis","Bacillus pumilus",
  "Bacillus megaterium","Bacillus circulans","Bacillus flexus",
  "Geobacillus stearothermophilus",
  "Pseudomonas fluorescens","Pseudomonas putida","Pseudomonas stutzeri",
  "Pseudomonas oleovorans","Pseudomonas veronii","Pseudomonas protegens",
  "Pseudomonas syringae","Pseudomonas alcaligenes","Pseudomonas mendocina",
  "Acinetobacter johnsonii","Acinetobacter lwoffii","Acinetobacter junii",
  "Acinetobacter radioresistens","Acinetobacter calcoaceticus",
  "Acinetobacter ursingii","Acinetobacter guillouiae",
  "Serratia marcescens","Serratia liquefaciens",
  # --- sequencing artefact -------------------------------------------------
  "Escherichia virus phiX174","Escherichia phage phiX174"
)

## ===========================================================================
## C. CONFLICTED GENERA
##    Contain BOTH documented contaminant species AND genuine gut residents.
##    A genus-level blacklist destroys real signal here; a species-level one
##    is only as good as the species call, which at <1x coverage is poor.
##    These are flagged "conflicted" and must be resolved at read level.
## ===========================================================================

CONFLICTED_GENERA <- c(
  "Escherichia","Shigella","Klebsiella","Enterobacter","Citrobacter",
  "Streptococcus","Staphylococcus","Pseudomonas","Acinetobacter","Bacillus",
  "Clostridium","Corynebacterium","Bifidobacterium","Lactobacillus",
  "Enterococcus","Serratia","Neisseria","Haemophilus","Actinomyces",
  "Rothia","Veillonella","Moraxella","Stenotrophomonas","Burkholderia",
  "Achromobacter","Micrococcus","Flavobacterium","Chryseobacterium"
)

## ===========================================================================
## D. GUT ALLOWLIST
##    Positive evidence of human colonic residency. This is the filter that
##    actually matters for near-sterile tissue: absence from a gut catalogue
##    is much stronger evidence of contamination than presence on a blacklist.
##
##    For production use, replace/augment with UHGG. Download the MGnify
##    human-gut catalogue metadata (genomes-all_metadata.tsv) from the MGnify
##    FTP under mgnify_genomes/human-gut/, take the species-level taxonomy
##    column, and pass the resulting character vector as `gut_catalogue`.
##    (Check the current path on the MGnify site; it moves between releases.)
## ===========================================================================

GUT_GENERA <- c(
  "Bacteroides","Phocaeicola","Parabacteroides","Alistipes","Prevotella",
  "Barnesiella","Odoribacter","Butyricimonas","Coprobacter","Paraprevotella",
  "Faecalibacterium","Roseburia","Agathobacter","Anaerobutyricum","Blautia",
  "Dorea","Coprococcus","Anaerostipes","Ruminococcus","Mediterraneibacter",
  "Lachnospira","Butyrivibrio","Subdoligranulum","Eubacterium","Oscillibacter",
  "Flavonifractor","Intestinibacter","Romboutsia","Turicibacter","Holdemanella",
  "Erysipelatoclostridium","Hungatella","Gemmiger","Ruthenibacterium",
  "Anaerotruncus","Fusicatenibacter","Marvinbryantia","Lachnoclostridium",
  "Bifidobacterium","Collinsella","Eggerthella","Gordonibacter","Adlercreutzia",
  "Slackia","Senegalimassilia","Olsenella",
  "Akkermansia","Victivallis","Methanobrevibacter","Methanosphaera",
  "Desulfovibrio","Bilophila","Sutterella","Parasutterella","Mailhella",
  "Fusobacterium","Leptotrichia","Campylobacter","Helicobacter",
  "Peptostreptococcus","Parvimonas","Anaerococcus","Finegoldia","Murdochiella",
  "Megasphaera","Dialister","Mitsuokella","Selenomonas","Acidaminococcus",
  "Phascolarctobacterium","Veillonella",
  "Christensenella","Catenibacterium","Solobacterium","Coprobacillus",
  "Bulleidia","Eisenbergiella","Emergencia","Anaerotignum","Terrisporobacter",
  "Sellimonas","Merdimonas","Faecalitalea","Clostridioides","Alloprevotella"
)

GUT_SPECIES <- c(
  # Bacteroidota
  "Bacteroides fragilis","Bacteroides thetaiotaomicron","Bacteroides uniformis",
  "Bacteroides ovatus","Bacteroides caccae","Bacteroides stercoris",
  "Bacteroides eggerthii","Bacteroides xylanisolvens","Bacteroides cellulosilyticus",
  "Bacteroides finegoldii","Bacteroides intestinalis","Bacteroides nordii",
  "Bacteroides salyersiae","Bacteroides clarus","Bacteroides fluxus",
  "Phocaeicola vulgatus","Bacteroides vulgatus","Phocaeicola dorei",
  "Bacteroides dorei","Phocaeicola massiliensis","Phocaeicola plebeius",
  "Parabacteroides distasonis","Parabacteroides merdae",
  "Parabacteroides goldsteinii","Parabacteroides johnsonii",
  "Alistipes putredinis","Alistipes shahii","Alistipes finegoldii",
  "Alistipes onderdonkii","Alistipes senegalensis","Alistipes indistinctus",
  "Alistipes ihumii","Prevotella copri","Segatella copri",
  "Barnesiella intestinihominis","Odoribacter splanchnicus",
  # Bacillota / Firmicutes
  "Faecalibacterium prausnitzii","Faecalibacterium duncaniae",
  "Roseburia intestinalis","Roseburia faecis","Roseburia hominis",
  "Roseburia inulinivorans",
  "Agathobacter rectalis","Eubacterium rectale","Eubacterium eligens",
  "Eubacterium ventriosum","Eubacterium siraeum","Eubacterium ramulus",
  "Anaerobutyricum hallii","Eubacterium hallii",
  "Blautia obeum","Blautia wexlerae","Blautia faecis","Blautia massiliensis",
  "Blautia hydrogenotrophica",
  "Dorea formicigenerans","Dorea longicatena",
  "Ruminococcus bromii","Ruminococcus gnavus","Mediterraneibacter gnavus",
  "Ruminococcus torques","Ruminococcus lactaris","Ruminococcus bicirculans",
  "Coprococcus comes","Coprococcus catus","Coprococcus eutactus",
  "Anaerostipes hadrus","Anaerostipes caccae",
  "Subdoligranulum variabile","Fusicatenibacter saccharivorans",
  "Flavonifractor plautii","Intestinibacter bartlettii",
  "Turicibacter sanguinis","Holdemanella biformis",
  "Erysipelatoclostridium ramosum","Hungatella hathewayi",
  "Clostridium bolteae","Clostridium symbiosum","Clostridium leptum",
  "Clostridium scindens","Clostridium innocuum","Clostridium perfringens",
  "Clostridium citroniae","Clostridioides difficile",
  "Enterococcus faecalis","Enterococcus faecium",
  "Lactobacillus rhamnosus","Lactobacillus ruminis","Limosilactobacillus reuteri",
  "Ligilactobacillus salivarius",
  "Megasphaera elsdenii","Dialister invisus","Dialister succinatiphilus",
  "Phascolarctobacterium faecium","Acidaminococcus intestini",
  "Catenibacterium mitsuokai",
  # Actinomycetota
  "Bifidobacterium longum","Bifidobacterium adolescentis","Bifidobacterium bifidum",
  "Bifidobacterium catenulatum","Bifidobacterium pseudocatenulatum",
  "Bifidobacterium breve","Bifidobacterium animalis",
  "Collinsella aerofaciens","Collinsella intestinalis",
  "Eggerthella lenta","Gordonibacter pamelaeae","Adlercreutzia equolifaciens",
  # Verrucomicrobiota / Euryarchaeota
  "Akkermansia muciniphila","Methanobrevibacter smithii",
  # Pseudomonadota
  "Escherichia coli","Klebsiella pneumoniae","Klebsiella oxytoca",
  "Klebsiella aerogenes","Enterobacter cloacae","Enterobacter hormaechei",
  "Citrobacter freundii","Citrobacter koseri","Proteus mirabilis",
  "Morganella morganii","Hafnia alvei",
  "Bilophila wadsworthia","Desulfovibrio piger","Sutterella wadsworthensis",
  "Parasutterella excrementihominis",
  # oral-derived but genuinely CRC-relevant in mucosa
  "Fusobacterium nucleatum","Fusobacterium periodonticum",
  "Peptostreptococcus stomatis","Peptostreptococcus anaerobius",
  "Parvimonas micra","Gemella morbillorum","Solobacterium moorei",
  "Streptococcus gallolyticus","Streptococcus anginosus","Streptococcus equinus",
  "Campylobacter concisus","Campylobacter ureolyticus","Campylobacter jejuni"
)

## ===========================================================================
## E. ORAL / UPPER-AIRWAY TAXA
##
## PROVENANCE OF THIS LIST -- read before citing it.
## The two lists below are a HAND-CURATED working set covering the commonly
## encountered oral genera and species. They were NOT extracted from a
## reference database, and they are not exhaustive (the oral microbiome has
## ~700 named taxa). For anything going into a paper, load the authoritative
## database instead, via load_oral_catalogue() below.
##
## Authoritative source:
##   expanded Human Oral Microbiome Database (eHOMD), https://www.homd.org
##   Escapa IF, Chen T, Huang Y, Gajare P, Dewhirst FE, Lemon KP (2018)
##     mSystems 3:e00187-18. doi:10.1128/mSystems.00187-18
##     770 species, 230 genera, curated, covering the whole aerodigestive tract.
##   Dewhirst FE et al. (2010) J Bacteriol 192:5002-5017. doi:10.1128/JB.00542-10
##     The founding ~600-taxon census behind HOMD.
##
## WHY "oral" IS SCORED NEUTRAL RATHER THAN PENALISED.
## Oral taxa in a colonic sample are not necessarily carryover. Schmidt TSB
## et al. (2019) eLife 8:e42693 (doi:10.7554/eLife.42693) found that oral-to-gut
## transmission is common in healthy people -- roughly one in three classifiable
## salivary cells also colonises the gut -- and is elevated in colorectal
## cancer. The opposing view is Rashidi A et al. (2021), who found no evidence
## of oral colonisation of the distal gut in healthy adults
## (PMID 34610963), so treat this as contested rather than settled.
##
## Several oral-origin taxa are among the most replicated CRC associations:
##   Wirbel J et al. (2019) Nat Med 25:679-689. doi:10.1038/s41591-019-0406-6
##     Meta-analysis of 8 studies; core set of 29 CRC-enriched species,
##     including Fusobacterium nucleatum, Peptostreptococcus stomatis and
##     Parvimonas micra. Also covers colibactin, fadA and bft directly.
## Penalising the whole oral category would therefore discard the findings
## most relevant to a colibactin study.
## ===========================================================================

ORAL_GENERA <- c(
  "Rothia","Neisseria","Haemophilus","Aggregatibacter","Capnocytophaga",
  "Porphyromonas","Tannerella","Treponema","Granulicatella","Abiotrophia",
  "Lautropia","Cardiobacterium","Eikenella","Kingella","Simonsiella",
  "Oribacterium","Stomatobaculum","Catonella","Johnsonella","Filifactor",
  "Peptidiphaga","Lachnoanaerobaculum","Alloscardovia","Scardovia",
  "Cryptobacterium","Mogibacterium","Atopobium","Schaalia"
)

ORAL_SPECIES <- c(
  "Streptococcus mitis","Streptococcus oralis","Streptococcus salivarius",
  "Streptococcus parasanguinis","Streptococcus sanguinis","Streptococcus gordonii",
  "Streptococcus cristatus","Streptococcus infantis","Streptococcus australis",
  "Streptococcus vestibularis","Streptococcus thermophilus",
  "Rothia mucilaginosa","Rothia dentocariosa","Rothia aeria",
  "Veillonella parvula","Veillonella atypica","Veillonella dispar",
  "Actinomyces odontolyticus","Actinomyces naeslundii","Actinomyces oris",
  "Prevotella melaninogenica","Prevotella intermedia","Prevotella nigrescens",
  "Porphyromonas gingivalis","Tannerella forsythia","Treponema denticola",
  "Haemophilus parainfluenzae","Neisseria flavescens","Neisseria subflava",
  "Leptotrichia buccalis","Leptotrichia wadei","Capnocytophaga gingivalis"
)

## ===========================================================================
## F. GENOTOXIN / COLIBACTIN-RELEVANT TAXA
##    Flagged separately because these are your a-priori hypotheses. The pks
##    (clb) island is not confined to E. coli -- it occurs in other
##    Enterobacteriaceae, so taxonomy alone will neither confirm nor exclude it.
## ===========================================================================

PKS_CAPABLE <- c(
  "Escherichia coli","Klebsiella pneumoniae","Klebsiella aerogenes",
  "Enterobacter aerogenes","Citrobacter koseri","Klebsiella oxytoca"
)

OTHER_GENOTOXIN <- c(
  "Bacteroides fragilis",        # fragilysin (ETBF)
  "Morganella morganii",         # indolimines
  "Escherichia albertii",        # CDT
  "Campylobacter jejuni",        # CDT
  "Helicobacter hepaticus",      # CDT
  "Shigella dysenteriae",        # CDT
  "Fusobacterium nucleatum",     # FadA / pro-tumourigenic
  "Enterococcus faecalis",       # extracellular superoxide
  "Peptostreptococcus anaerobius","Parvimonas micra"
)

## ===========================================================================
## G. FUNCTIONS
## ===========================================================================

#' Normalise a MetaPhlAn taxon token to a plain name.
#' "s__Escherichia_coli" -> "Escherichia coli"
mp_clean <- function(x) {
  x <- sub("^[a-z]__", "", x)
  x <- gsub("_", " ", x)
  # NCBI brackets a genus whose placement is disputed: "[Ruminococcus] gnavus",
  # "[Clostridium] scindens", "[Eubacterium] hallii". These are REAL and often
  # important gut organisms; leaving the brackets in makes the genus
  # "[Ruminococcus]", which matches nothing. Strip them.
  x <- gsub("[][]", "", x)
  trimws(gsub("  +", " ", x))
}

#' Pull the token for one rank out of lineage strings.
#' Returns NA where that rank is absent.
mp_rank_token <- function(lineage, rank) {
  pat <- paste0("(^|\\|)", rank, "__([^|]*)")
  m <- regexpr(pat, lineage)
  out <- rep(NA_character_, length(lineage))
  hit <- m > 0
  # regmatches() returns only the elements that matched, so assign into out[hit]
  out[hit] <- sub(paste0("^\\|?", rank, "__"), "", regmatches(lineage, m))
  out
}

#' Build a proper taxon name from a lineage string.
#'
#' Some profilers write the full binomial in the species field
#' ("s__Halopenitus_persicus"); others write only the epithet
#' ("s__persicus") and leave the genus in the g__ field. The second form is
#' unmatchable against any species list, and bare epithets are not even unique
#' (s__marina occurs under several genera), so this reconstructs
#' "Halopenitus persicus" from "g__Halopenitus|s__persicus".
#'
#' @param lineage character vector of full lineage strings
#' @param rank    the rank to name ("s", "g", "t", ...)
mp_name_from_lineage <- function(lineage, rank = "s") {
  tgt <- mp_rank_token(lineage, rank)
  nm  <- mp_clean(ifelse(is.na(tgt), NA_character_, tgt))
  if (rank != "s") return(nm)
  gen <- mp_clean(mp_rank_token(lineage, "g"))
  need_g <- !is.na(nm) & !is.na(gen) & nm != gen &
    !startsWith(nm, paste0(gen, " "))
  nm[need_g] <- paste(gen[need_g], nm[need_g])
  nm
}

#' Build a taxon x sample count matrix with correct binomial rownames.
#'
#' Repairs genus-dropped species fields and sums any rows that collapse to the
#' same name once repaired.
#'
#' @param counts  numeric matrix or data.frame, rows aligned to `lineage`
#' @param lineage character vector of full lineage strings
#' @param rank    single-letter rank the rows were already filtered to
#' @return list(mat = integer matrix, lineage = character)
build_taxon_matrix <- function(counts, lineage, rank = "s") {
  m <- as.matrix(counts)
  m[is.na(m)] <- 0
  m <- round(m)
  nm <- mp_name_from_lineage(lineage, rank)
  bad <- is.na(nm) | !nzchar(nm)
  if (any(bad)) nm[bad] <- mp_clean(sub(".*\\|", "", lineage[bad]))
  if (anyDuplicated(nm)) {
    message("build_taxon_matrix(): summing ", sum(duplicated(nm)),
            " row(s) collapsing to a duplicate taxon name.")
    first   <- !duplicated(nm)
    m       <- rowsum(m, group = nm, reorder = FALSE)
    lineage <- lineage[first]
    nm      <- nm[first]
  }
  rownames(m) <- nm
  list(mat = m, lineage = lineage)
}

#' Is this an unnamed MetaPhlAn 4 SGB/GGB placeholder?
#' These have no literature, no catalogue entry and no interpretation at low
#' depth. Treat them as their own category, not as "unknown but plausible".
is_unnamed <- function(x) {
  grepl("^(GGB|SGB|FGB|OFGB|CFGB|PFGB|KFGB)[0-9]", x) |
    grepl("(GGB|SGB)[0-9]+", x) & !grepl("^[A-Z][a-z]+ [a-z]", x)
}

#' Annotate a vector of taxon names.
#'
#' @param taxa          character vector, either MetaPhlAn tokens
#'                      ("s__Escherichia_coli") or plain names.
#' @param rank          "s" (species) or "g" (genus).
#' @param gut_catalogue optional character vector of species names from an
#'                      external catalogue (e.g. UHGG). If supplied, it takes
#'                      precedence over the built-in GUT_SPECIES list.
#' @param lineage       optional character vector of FULL lineage strings
#'                      ("d__Bacteria|p__...|s__..."), same length and order as
#'                      `taxa`. Strongly recommended: it lets whole
#'                      environmental clades be excluded without enumerating
#'                      every genus, which is what a full-NCBI database needs.
#' @return data.frame, one row per input taxon, in input order.
annotate_taxa <- function(taxa, rank = c("s", "g"), gut_catalogue = NULL,
                          lineage = NULL, oral_catalogue = NULL) {
  rank <- match.arg(rank)
  nm   <- mp_clean(taxa)
  
  ## Safety net for tables that drop the genus from the species field
  ## ("s__persicus" rather than "s__Halopenitus_persicus"). A bare epithet is
  ## lower-case and has no space; if a lineage is available, rebuild the
  ## binomial from the g__ token. Without the lineage nothing can be done, so
  ## warn rather than silently return an all-"clean" table.
  if (rank == "s") {
    bare <- !is.na(nm) & nzchar(nm) & !grepl("^[A-Z]", nm) & !grepl(" ", nm)
    if (any(bare)) {
      if (!is.null(lineage)) {
        fixed <- mp_name_from_lineage(lineage, "s")
        nm[bare & !is.na(fixed)] <- fixed[bare & !is.na(fixed)]
        message("annotate_taxa(): rebuilt ", sum(bare),
                " genus-dropped species name(s) from the lineage.")
      } else {
        warning("annotate_taxa(): ", sum(bare), " taxa look like bare species ",
                "epithets (e.g. '", nm[which(bare)[1]], "'). Pass `lineage=` ",
                "so the genus can be recovered, or nothing will match.",
                call. = FALSE)
      }
    }
  }
  
  unnamed <- is_unnamed(nm)
  
  # NB: do NOT use ifelse() here. `rank` is length 1, and ifelse() returns a
  # result the same length as its TEST, so ifelse(rank=="s", vec1, vec2) would
  # silently return a single element -- which then recycles through every
  # %in% check and makes every taxon look "clean".
  if (rank == "s") {
    genus   <- sub(" .*", "", nm)
    species <- nm
  } else {
    genus   <- nm
    species <- rep(NA_character_, length(nm))
  }
  
  # "Rhizobium sp. ACO-34A", "Pseudomonas sp. WL3", "uncultured bacterium":
  # NCBI/Bracken strain-level entries with no formal species description.
  undescribed <- grepl(" sp\\.| sp$|uncultured|unclassified| bacterium( |$)",
                       nm, ignore.case = TRUE)
  
  gut_sp <- if (!is.null(gut_catalogue)) unique(c(GUT_SPECIES, gut_catalogue)) else GUT_SPECIES
  
  # eHOMD, if supplied. Species are taken wholesale; GENERA are added only
  # where they are not already gut or contaminant genera, because ~230 eHOMD
  # genera include many (Streptococcus, Prevotella, Veillonella, Fusobacterium)
  # that are perfectly at home in the colon. Applying those at genus level
  # would relabel real gut organisms as oral.
  oral_sp <- ORAL_SPECIES; oral_gn <- ORAL_GENERA
  if (!is.null(oral_catalogue)) {
    oral_sp <- unique(c(oral_sp, oral_catalogue$species))
    oral_gn <- unique(c(oral_gn, setdiff(oral_catalogue$genera,
                                         c(GUT_GENERA, GENUS_CONTAM))))
  }
  
  ## --- whole-clade environmental exclusion (needs the lineage) -------------
  env_clade <- rep(FALSE, length(nm))
  if (!is.null(lineage)) {
    if (length(lineage) != length(taxa))
      stop("`lineage` must be the same length as `taxa`.")
    env_clade <- grepl(paste(ENV_CLADE_PATTERNS, collapse = "|"), lineage)
  }
  
  ## --- blacklist ------------------------------------------------------------
  bl <- rep("clean", length(nm))
  if (rank == "s") {
    bl[genus %in% GENUS_CONTAM]     <- "contaminant"
    bl[species %in% SPECIES_CONTAM] <- "contaminant"
    # At SPECIES rank a conflicted genus is only genuinely conflicted when the
    # species itself has gut evidence. An obscure environmental congener
    # (Pseudomonas veronii, say) is not "ambiguous" -- it simply has no gut
    # support, and the allowlist already excludes it. Flagging every member of
    # a big genus as conflicted made 54% of a full-NCBI profile unreadable.
    bl[genus %in% CONFLICTED_GENERA & species %in% gut_sp] <- "conflicted_gut_plausible"
  } else {
    bl[genus %in% GENUS_CONTAM]      <- "contaminant"
    bl[genus %in% CONFLICTED_GENERA] <- "conflicted"
  }
  bl[env_clade] <- "contaminant"
  bl[unnamed]   <- "unknown"
  
  ## --- allowlist ------------------------------------------------------------
  al <- rep("unknown", length(nm))
  if (rank == "s") {
    al[genus   %in% oral_gn] <- "oral"
    al[species %in% oral_sp] <- "oral"
    al[genus   %in% GUT_GENERA] <- "gut"
    al[species %in% gut_sp]     <- "gut"
  } else {
    al[genus %in% oral_gn]     <- "oral"
    al[genus %in% GUT_GENERA]  <- "gut"
  }
  # "environmental" is more informative than "unknown": it means we have
  # positive evidence this CANNOT be a colonic resident, not merely absence
  # of evidence that it is.
  al[env_clade | genus %in% GENUS_CONTAM] <- "environmental"
  al[unnamed] <- "unnamed_sgb"
  
  ## --- a-priori hypothesis flags -------------------------------------------
  pks <- if (rank == "s") species %in% PKS_CAPABLE else
    genus %in% unique(sub(" .*", "", PKS_CAPABLE))
  gtx <- if (rank == "s") species %in% OTHER_GENOTOXIN else
    genus %in% unique(sub(" .*", "", OTHER_GENOTOXIN))
  
  out <- data.frame(
    taxon             = taxa,
    name              = nm,
    genus             = genus,
    species           = species,
    unnamed_sgb       = unnamed,
    undescribed_strain = undescribed,
    env_clade         = env_clade,
    blacklist_status  = bl,
    allowlist_status  = al,
    pks_capable       = pks,
    other_genotoxin   = gtx,
    stringsAsFactors  = FALSE
  )
  
  # Loud failure instead of a silent all-"clean" table. If a real metagenome
  # matches nothing, the lists are not being applied and the result is garbage.
  if (length(nm) > 50 && all(bl == "clean"))
    warning("annotate_taxa(): NOTHING matched any list across ", length(nm),
            " taxa. Check that the names look like '", nm[1],
            "' -- if they carry a rank prefix or odd separators, mp_clean() ",
            "is not stripping them. Run annotation_report() to inspect.",
            call. = FALSE)
  out
}

#' Print what actually matched, and the biggest taxa that matched nothing.
#' Run this every time -- it is how you notice a list that isn't firing.
annotation_report <- function(ann, mat = NULL, top_n = 30) {
  cat("\n--- annotation summary (", nrow(ann), " taxa ) ---\n", sep = "")
  print(table(ann$blacklist_status))
  print(table(ann$allowlist_status))
  cat("pks-capable:", sum(ann$pks_capable),
      " other genotoxin:", sum(ann$other_genotoxin),
      " undescribed 'sp.' strains:", sum(ann$undescribed_strain), "\n")
  if (!is.null(mat)) {
    stopifnot(nrow(ann) == nrow(mat))
    tot <- rowSums(mat)
    unk <- which(ann$allowlist_status == "unknown" & ann$blacklist_status == "clean")
    if (length(unk)) {
      unk <- unk[order(-tot[unk])][seq_len(min(top_n, length(unk)))]
      cat("\nMost abundant taxa matching NEITHER list (",
          round(100 * sum(tot[ann$allowlist_status == "unknown" &
                                ann$blacklist_status == "clean"]) / sum(tot), 1),
          "% of reads).\nIf these dominate, add their genera to the lists in",
          " taxon_annotation.R:\n", sep = "")
      print(data.frame(taxon = ann$name[unk], genus = ann$genus[unk],
                       reads = tot[unk], row.names = NULL))
    }
  }
  invisible(NULL)
}

#' Add data-driven contamination evidence from your own count matrix.
#'
#' @param ann      output of annotate_taxa(), rows aligned to rownames(mat)
#' @param mat      taxa x samples integer count matrix
#' @param blanks   optional character vector of column names that are true
#'                 negative controls. If given, blank-based evidence dominates.
#' @param genome_sizes named numeric vector (bp), optional
score_evidence <- function(ann, mat, blanks = NULL, read_len = 150,
                           paired = TRUE, genome_sizes = NULL,
                           min_reads = 10, default_genome = 4e6) {
  stopifnot(nrow(ann) == nrow(mat))
  is_blank <- if (!is.null(blanks)) colnames(mat) %in% blanks else rep(FALSE, ncol(mat))
  smat     <- mat[, !is_blank, drop = FALSE]
  libsize  <- colSums(smat)
  
  gs <- rep(default_genome, nrow(mat))
  if (!is.null(genome_sizes)) {
    # match on full name first, fall back to genus. Indexing a named vector
    # with an absent name yields NA (not NULL), so coalesce explicitly.
    hit <- unname(genome_sizes[ann$name])
    alt <- unname(genome_sizes[ann$genus])
    hit[is.na(hit)] <- alt[is.na(hit)]
    gs[!is.na(hit)] <- hit[!is.na(hit)]
  }
  bases <- read_len * (if (paired) 2 else 1)
  max_cov <- apply(smat * bases / gs, 1, max)
  
  relab <- sweep(smat, 2, pmax(libsize, 1), "/")
  rho <- apply(relab, 1, function(x) {
    if (sum(x > 0) < 4) return(NA_real_)
    suppressWarnings(stats::cor(log10(libsize + 1), x, method = "spearman"))
  })
  
  prev <- rowMeans(smat >= min_reads)
  
  blank_prev <- if (any(is_blank))
    rowMeans(mat[, is_blank, drop = FALSE] >= min_reads) else rep(NA_real_, nrow(mat))
  # ratio of mean relative abundance in samples vs blanks; <1 means the blanks
  # have MORE of it than the samples do
  blank_ratio <- if (any(is_blank)) {
    bm <- mat[, is_blank, drop = FALSE]
    bl_rel <- rowMeans(sweep(bm, 2, pmax(colSums(bm), 1), "/"))
    signif(rowMeans(relab) / pmax(bl_rel, 1e-9), 3)
  } else rep(NA_real_, nrow(mat))
  
  ## --- additive evidence score ---------------------------------------------
  sc <- rep(0, nrow(mat))
  sc <- sc + ifelse(ann$allowlist_status == "gut", 2, 0)
  # "oral" is deliberately scored 0, not negative: F. nucleatum, P. micra and
  # P. stomatis are oral in origin but are real, replicated CRC mucosal findings.
  sc <- sc - ifelse(ann$allowlist_status == "unnamed_sgb", 1, 0)
  sc <- sc - ifelse(ann$allowlist_status == "environmental", 2, 0)
  sc <- sc - ifelse(ann$blacklist_status == "contaminant", 3, 0)
  sc <- sc + ifelse(max_cov >= 0.05, 1, 0)
  sc <- sc + ifelse(max_cov >= 0.50, 1, 0)
  sc <- sc + ifelse(prev >= 0.10, 1, 0)
  sc <- sc + ifelse(!is.na(rho) & rho >= -0.1, 1, 0)
  sc <- sc - ifelse(!is.na(rho) & rho < -0.3, 2, 0)
  if (any(is_blank)) {
    sc <- sc - ifelse(blank_prev > 0.5, 3, 0)
    sc <- sc + ifelse(blank_prev == 0, 2, 0)
    sc <- sc + ifelse(!is.na(blank_ratio) & blank_ratio > 10, 1, 0)
  }
  
  tier <- cut(sc, breaks = c(-Inf, -1, 1, 3, Inf),
              labels = c("D_discard", "C_weak", "B_plausible", "A_confident"))
  
  out <- cbind(ann, data.frame(
    max_coverage_X   = signif(max_cov, 3),
    prevalence       = round(prev, 3),
    rho_vs_depth     = round(rho, 2),
    blank_prevalence = round(blank_prev, 3),
    sample_blank_ratio = blank_ratio,
    evidence_score   = sc,
    tier             = as.character(tier),
    stringsAsFactors = FALSE
  ))
  # single convenience flag, but keep the components visible
  out$USE_IN_ANALYSIS <- out$tier %in% c("A_confident", "B_plausible")
  out
}

#' Load the eHOMD oral taxon list from a downloaded HOMD table.
#'
#' Replaces the hand-curated ORAL_SPECIES/ORAL_GENERA with the curated
#' database, so the oral category becomes citable rather than ad hoc.
#'
#' Get the file from https://www.homd.org (Downloads -> taxonomy table). Any
#' delimited table with columns giving genus and species will work; the
#' function looks for the usual HOMD column names and falls back to guessing.
#'
#' @return list(species = <character>, genera = <character>), or NULL.
load_oral_catalogue <- function(path) {
  if (!file.exists(path)) {
    warning("Oral catalogue not found at ", path, " -- using built-in list.")
    return(NULL)
  }
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  d <- utils::read.delim(path, sep = sep, header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE)
  gcol <- intersect(c("Genus", "genus"), names(d))
  scol <- intersect(c("Species", "species"), names(d))
  full <- intersect(c("Name", "name", "Taxon", "scientific_name"), names(d))
  sp <- if (length(gcol) && length(scol)) {
    paste(trimws(d[[gcol[1]]]), trimws(d[[scol[1]]]))
  } else if (length(full)) {
    trimws(d[[full[1]]])
  } else {
    stop("Could not find genus/species or name columns in ", path)
  }
  sp <- sp[nzchar(sp) & !grepl("\\bsp\\.?$|unclassified|uncultured", sp,
                               ignore.case = TRUE)]
  list(species = unique(sp), genera = unique(sub(" .*", "", sp)))
}

#' Convenience: load an external gut catalogue species list from a file.
#' Accepts UHGG genomes-all_metadata.tsv (uses the Lineage column) or a plain
#' one-species-per-line text file.
load_gut_catalogue <- function(path) {
  if (!file.exists(path)) {
    warning("Gut catalogue not found at ", path, " -- using built-in list only.")
    return(NULL)
  }
  first <- readLines(path, n = 1)
  if (grepl("\t", first)) {
    d <- utils::read.delim(path, sep = "\t", header = TRUE,
                           check.names = FALSE, stringsAsFactors = FALSE)
    col <- intersect(c("Lineage", "lineage", "GTDB_taxonomy", "taxonomy"), names(d))
    if (!length(col)) stop("Could not find a lineage/taxonomy column in ", path)
    lin <- d[[col[1]]]
    sp  <- sub(".*;s__", "", lin)
    sp  <- sp[nzchar(sp) & !grepl(";", sp)]
    unique(trimws(sp))
  } else {
    unique(trimws(readLines(path)))
  }
}

## ===========================================================================
## H. SELF-TEST
##    Run test_annotation() after source()-ing this file. It asserts that each
##    list actually fires. Every class of bug that produced an all-"clean"
##    table would fail here.
## ===========================================================================

test_annotation <- function() {
  probe <- c("s__Escherichia_coli", "s__Ralstonia_pickettii",
             "s__Rhizobium_sp._ACO-34A", "s__Faecalibacterium_prausnitzii",
             "s__Cutibacterium_acnes", "s__Bacteroides_fragilis",
             "s__GGB9710_SGB15265", "s__Fusobacterium_nucleatum",
             "s__Staphylococcus_epidermidis", "s__Klebsiella_pneumoniae",
             "s__[Ruminococcus]_gnavus", "s__[Clostridium]_scindens",
             "s__Stutzerimonas_stutzeri", "s__Natronomonas_moolapensis",
             "s__Pseudomonas_aeruginosa")
  lin <- c(rep("d__Bacteria|p__x", 13),
           "d__Archaea|p__Halobacteriota|c__Halobacteria|s__Natronomonas_moolapensis",
           "d__Bacteria|p__x")
  a <- annotate_taxa(probe, rank = "s", lineage = lin)
  print(a[, c("name", "genus", "blacklist_status", "allowlist_status",
              "pks_capable", "undescribed_strain")])
  g <- function(n, col) a[[col]][match(n, a$name)]
  stopifnot(
    nrow(a) == length(probe),
    length(unique(a$genus)) > 1,   # the all-"clean" bug collapsed this to 1
    g("Ralstonia pickettii",  "blacklist_status") == "contaminant",
    g("Cutibacterium acnes",  "blacklist_status") == "contaminant",
    g("Rhizobium sp. ACO-34A","blacklist_status") == "contaminant",
    g("Rhizobium sp. ACO-34A","undescribed_strain"),
    g("Staphylococcus epidermidis", "blacklist_status") == "contaminant",
    g("Stutzerimonas stutzeri", "blacklist_status") == "contaminant",
    # excluded by LINEAGE (c__Halobacteria), not by an enumerated genus
    g("Natronomonas moolapensis", "env_clade"),
    g("Natronomonas moolapensis", "allowlist_status") == "environmental",
    # NCBI brackets must be stripped or these real gut taxa are invisible
    g("Ruminococcus gnavus",  "allowlist_status") == "gut",
    g("Clostridium scindens", "allowlist_status") == "gut",
    # a conflicted genus is only "conflicted" when the SPECIES has gut evidence
    g("Escherichia coli",     "blacklist_status") == "conflicted_gut_plausible",
    g("Pseudomonas aeruginosa", "blacklist_status") == "clean",
    g("Escherichia coli",     "allowlist_status") == "gut",
    g("Escherichia coli",     "pks_capable"),
    g("Klebsiella pneumoniae","pks_capable"),
    g("Faecalibacterium prausnitzii", "allowlist_status") == "gut",
    g("Bacteroides fragilis", "other_genotoxin"),
    g("GGB9710 SGB15265",     "allowlist_status") == "unnamed_sgb",
    g("Fusobacterium nucleatum", "allowlist_status") == "gut"
  )
  ## --- genus-dropped species fields (GTDB-style lineages) -----------------
  ## Some profiler versions write only the epithet in s__, leaving the genus
  ## in g__. Bare epithets match nothing and are not unique across genera.
  lin2 <- c(
    "d__Archaea|k__Methanobacteriati|p__Methanobacteriota|c__Halobacteria|o__Halobacteriales|f__Haloferacaceae|g__Halopenitus|s__persicus",
    "d__Archaea|k__Methanobacteriati|p__Methanobacteriota|c__Halobacteria|o__Halobacteriales|f__Natronomonadaceae|g__Natronomonas|s__moolapensis",
    "d__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__coli",
    "d__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Faecalibacterium|s__prausnitzii",
    # already-full binomial must NOT be doubled up
    "d__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Roseburia|s__Roseburia_intestinalis"
  )
  nm2 <- mp_name_from_lineage(lin2, "s")
  print(data.frame(lineage_tail = sub(".*\\|", "", lin2), rebuilt = nm2))
  stopifnot(
    nm2[1] == "Halopenitus persicus",
    nm2[2] == "Natronomonas moolapensis",
    nm2[3] == "Escherichia coli",
    nm2[4] == "Faecalibacterium prausnitzii",
    nm2[5] == "Roseburia intestinalis"            # not doubled
  )
  b <- annotate_taxa(sub(".*\\|", "", lin2), rank = "s", lineage = lin2)
  stopifnot(
    b$allowlist_status[1] == "environmental",     # c__Halobacteria
    b$allowlist_status[2] == "environmental",
    b$env_clade[1], b$env_clade[2],
    !b$env_clade[3],                              # p__Methanobacteriota is NOT
    b$allowlist_status[3] == "gut",               # blanket-excluded, since
    b$allowlist_status[4] == "gut",               # M. smithii lives there
    b$pks_capable[3]
  )
  
  # bare epithets must collapse to the same name only via the genus
  cm <- build_taxon_matrix(matrix(1:10, nrow = 5), lin2, "s")
  stopifnot(nrow(cm$mat) == 5, !anyDuplicated(rownames(cm$mat)))
  
  cat("\ntest_annotation(): PASSED -- all lists are firing.\n")
  invisible(a)
}

message("taxon_annotation.R loaded: ",
        length(SPECIES_CONTAM), " contaminant species, ",
        length(GENUS_CONTAM), " contaminant genera, ",
        length(GUT_SPECIES), " gut species, ",
        length(CONFLICTED_GENERA), " conflicted genera. ",
        "Run test_annotation() to verify.")