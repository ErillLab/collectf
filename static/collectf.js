$(document).ready(function() {
// onclick
    if($("#id_1-TF_species_same").length > 0) {
	$("#id_1-TF_species_same").click(function(){
	    $("#id_1-TF_species")[0].disabled = this.checked;
	});

	$("#id_1-site_species_same").click(function(){
	    $("#id_1-site_species")[0].disabled = this.checked;
	});

	// when page is loaded, if checkboxes are checked
	$("#id_1-TF_species")[0].disabled = $("#id_1-TF_species_same")[0].checked;
	$("#id_1-site_species")[0].disabled = $("#id_1-site_species_same")[0].checked;
    }

    // twitter bootstrap popover
    $('body').popover({
	selector: '[data-toggle="popover"]',
	trigger: 'hover',
	placement: 'right',
	html: 'true',
    });

    $('body').tooltip({
	selector: '[data-toggle="tooltip"]',
	trigger: 'hover',
	placement: 'top',
	html: 'true',
    });
    
    

    // scroll page down a little bit
    //window.addEventListener("hashchange", function() { scrollBy(0, -50); });

    defaultSpeed = 200;

    // JS for reported sites step of curation form.
    if ($('#id_3-is_chip_data').not(':checked')) {
	$('#id_3-assay_conditions').parent().parent().hide(defaultSpeed);
	$('#id_3-chip_method_notes').parent().parent().hide(defaultSpeed);
	$('#id_3-chip_data_extra_field').parent().parent().hide(defaultSpeed);
    }
    if ($('#id_3-is_chip_data').is(':checked')) {
	$('#id_3-assay_conditions').parent().parent().show(defaultSpeed);
	$('#id_3-chip_method_notes').parent().parent().show(defaultSpeed);
	$('#id_3-chip_data_extra_field').parent().parent().show(defaultSpeed);
    }
    $("#id_3-is_chip_data").on('change', function() {
	if ($(this).is(":checked")) {
	    $("#id_3-assay_conditions").parent().parent().show(defaultSpeed);
	    $("#id_3-chip_method_notes").parent().parent().show(defaultSpeed);
	    $('#id_3-chip_data_extra_field').parent().parent().show(defaultSpeed);
	} else {
	    $("#id_3-assay_conditions").parent().parent().hide(defaultSpeed);
	    $("#id_3-chip_method_notes").parent().parent().hide(defaultSpeed);
	    $('#id_3-chip_data_extra_field').parent().parent().hide(defaultSpeed);
	}
    });

    if ($('#id_3-has_quantitative_data').not(':checked')) {
	$('#id_3-quantitative_data_format').parent().parent().hide(defaultSpeed);
    }
    if ($('#id_3-has_quantitative_data').is(':checked')) {
	$('#id_3-quantitative_data_format').parent().parent().show(defaultSpeed);
    }
    $('#id_3-has_quantitative_data').on('change', function() {
	if ($(this).is(':checked')) {
	    $('#id_3-quantitative_data_format').parent().parent().show(defaultSpeed);
	} else {
	    $('#id_3-quantitative_data_format').parent().parent().hide(defaultSpeed);
	}
    });

});


