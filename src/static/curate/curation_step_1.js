// Javascript utilities for Curation step 1: Genome and TF accesion numbers

function disableTFSpecies() {
    // Disables TF-species field if it is checked as same.
    var cb = $("#id_1-TF_species_same");
    var text = $("#id_1-TF_species");
    text.attr("disabled", cb.attr("checked"));
    cb.click(function() {text.attr("disabled", this.checked);});
}

function disableSiteSpecies() {
    // Disables site-species field if it is checked as same.
    var cb = $("#id_1-site_species_same");
    var text = $("#id_1-site_species");
    text.attr("disabled", cb.attr("checked"));
    cb.click(function() {text.attr("disabled", this.checked);});
}

function typeaheadUpdater(item) {
    // Function to be called when typeahead input element is selected.
    return item.split(' ')[0];
}

function genomeAccessionTypeAhead() {
    // Typeahead for genome accession number textbox.
    genomeFields = $("[id^=id_1-genome_accession]");
    genomeFields.attr('autocomplete', 'off');
    $.getJSON("/curate/json_get_genomes", function (data) {
        var genomes = $.map(data, function (n, i) {
            return n['genome_accession'] + ' - ' + n['organism'];
        });
        genomeFields.typeahead({
            source: genomes,
            updater: typeaheadUpdater
        });
    });
}
function TFAccesssionTypeAhead() {
    // Typeahead for TF accession number textbox
    TFFields = $("[id^=id_1-TF_accession]");
    TFFields.attr('autocomplete', 'off');
    $.getJSON("/curate/json_get_TF_instances", function (data) {
        var TFs = $.map(data, function (n, i) {
            return (n['uniprot_accession'] + ' (' + n['refseq_accession'] +
                    ') ' + n['description']);
        });
        TFFields.typeahead({
            source: TFs,
            updater: typeaheadUpdater
        });
    });
}

function addToggleLinks() {
    // Adds "Extra genome accession" and "Extra TF accession" toggle links to
    // the help text.
    var s = $("#id_1-genome_accession").next().html();
    var add_text = ' <small><a href="#" onclick="toggleGenomeFields(); return false;">Toggle extra genome accession fields</a></small>';
    $("#id_1-genome_accession").next().html(s + add_text);
    // Same thing for TF accession fields
    var s = $("#id_1-TF_accession").next().html();
    var add_text = ' <small><a href="#" onclick="toggleTFFields(); return false;">Toggle extra TF accession fields</a></small>';
    $("#id_1-TF_accession").next().html(s + add_text);
}

function toggleGenomeFields() {
    var genome_accessions = $("[id^='id_1-genome_accession_']");
    for (i=0; i<genome_accessions.length; i++) {
        var g = $(genome_accessions[i]);
        if(!(g.val())) {
            g.closest('.form-group').toggle('fast');
        }
    }
}

function toggleTFFields() {
    // hide extra genome accession fields
    var TF_accessions = $("[id^='id_1-TF_accession_']");
    for (i=0; i<TF_accessions.length; i++) {
        var t = $(TF_accessions[i]);
        if (!(t.val())) {
            t.closest('.form-group').toggle('fast');
        }
    }
}

$(document).ready(function() {
    disableTFSpecies();
    disableSiteSpecies();
    genomeAccessionTypeAhead();
    TFAccesssionTypeAhead();
    addToggleLinks();
    toggleGenomeFields();
    toggleTFFields();
});
