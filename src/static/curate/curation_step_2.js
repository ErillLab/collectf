function externalDBToggle() {
    // Collapses external DB type and accession number fields together
    var doms =  $("[id^='id_2-external_db_accession']");
    for(var i=0; i<doms.length; i++) {
        var d = $(doms[i]);
        console.log(d.val());
        if (!(d.val())) {
            d.parent().parent().toggle('fast');
            d.parent().parent().prev().toggle('fast');
        }
    }
}

function complexNotesToggle() {
    // Collapses complex notes
    var dom_2 = $('#id_2-complex_notes');
    if (!dom_2.val()) {
        dom_2.parent().parent().toggle('fast');
    }
}

$(document).ready(function() {
    externalDBToggle();
    $('#id_2-has_external_db').change(externalDBToggle);
    complexNotesToggle();
    $('#id_2-forms_complex').change(complexNotesToggle);
    // Change the ul style (hide bullet points on the list of experimental techniques)
    $('ul').css('list-style-type', 'none');
    // Add title "Additional information" below "Experimental process" field
    $('.form-group:has(#id_2-forms_complex)').before("<div class='page-header'><h3>Additional information</h3></div>");
});

