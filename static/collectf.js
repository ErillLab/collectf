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
	$('#id_3-peak_calling_method').parent().parent().hide(defaultSpeed);
	$('#id_3-assay_conditions').parent().parent().hide(defaultSpeed);
	$('#id_3-chip_method_notes').parent().parent().hide(defaultSpeed);
	$('#id_3-chip_data_extra_field').parent().parent().hide(defaultSpeed);
    }
    if ($('#id_3-is_chip_data').is(':checked')) {
	$('#id_3-peak_calling_method').parent().parent().show(defaultSpeed);
	$('#id_3-assay_conditions').parent().parent().show(defaultSpeed);
	$('#id_3-chip_method_notes').parent().parent().show(defaultSpeed);
	$('#id_3-chip_data_extra_field').parent().parent().show(defaultSpeed);
    }
    $("#id_3-is_chip_data").on('change', function() {
	if ($(this).is(":checked")) {
	    $("#id_3-peak_calling_method").parent().parent().show(defaultSpeed);
	    $("#id_3-assay_conditions").parent().parent().show(defaultSpeed);
	    $("#id_3-chip_method_notes").parent().parent().show(defaultSpeed);
	    $('#id_3-chip_data_extra_field').parent().parent().show(defaultSpeed);
	} else {
	    $("#id_3-peak_calling_method").parent().parent().hide(defaultSpeed);
	    $("#id_3-assay_conditions").parent().parent().hide(defaultSpeed);
	    $("#id_3-chip_method_notes").parent().parent().hide(defaultSpeed);
	    $('#id_3-chip_data_extra_field').parent().parent().hide(defaultSpeed);
	}
    });

    /*
    // Hide <is_it_sequence> field by default
    if ($("#id_3-motif_associated").is(":checked")) {
	$("#id_3-is_coordinate").parent().parent().hide(defaultSpeed);
	$("#id_3-is_coordinate").prop('checked', false);
	$("#id_3-is_chip_seq_data").parent().parent().hide(defaultSpeed);
	$("#id_3-is_chip_seq_data").prop('checked', false);
	$("#id_3-peak_calling_method").parent().parent().hide(defaultSpeed);
	$("#id_3-assay_conditions").parent().parent().hide(defaultSpeed);
	$("#id_3-method_notes").parent().parent().hide(defaultSpeed);
    }

    
    if ($("#id_3-is_chip_seq_data").not(":checked")) {
	$("#id_3-peak_intensity").parent().parent().hide(defaultSpeed);
	$("#id_3-peak_calling_method").parent().parent().hide(defaultSpeed);
	$("#id_3-assay_conditions").parent().parent().hide(defaultSpeed);
	$("#id_3-method_notes").parent().parent().hide(defaultSpeed);
    }

    if ($("#id_3-is_chip_seq_data").is(":checked")) {
	$("#id_3-peak_intensity").parent().parent().show(defaultSpeed);
	$("#id_3-peak_calling_method").parent().parent().show(defaultSpeed);
	$("#id_3-assay_conditions").parent().parent().show(defaultSpeed);
	$("#id_3-method_notes").parent().parent().show(defaultSpeed);
    }

    

    // Based on whether <is_motif_associated> field is checked or not, hide/show
    // <is_it_sequence> field.
    $("#id_3-motif_associated").on('change', function() {
	if ($(this).is(":checked")) {
	    $("#id_3-is_coordinate").parent().parent().hide(defaultSpeed);
	    $("#id_3-is_chip_seq_data").parent().parent().hide(defaultSpeed);
	    $("#id_3-peak_calling_method").parent().parent().hide(defaultSpeed);
	    $("#id_3-assay_conditions").parent().parent().hide(defaultSpeed);
	    $("#id_3-method_notes").parent().parent().hide(defaultSpeed);
	} else {
	    $("#id_3-is_coordinate").parent().parent().show(defaultSpeed);
	    $("#id_3-is_chip_seq_data").parent().parent().show(defaultSpeed);

	}
    });

    $("#id_3-is_chip_seq_data").on('change', function() {
	if ($(this).is(":checked")) {
	    $("#id_3-is_coordinate").prop("checked", true);
	    $("#id_3-peak_intensity").parent().parent().show(defaultSpeed);
	    $("#id_3-peak_calling_method").parent().parent().show(defaultSpeed);
	    $("#id_3-assay_conditions").parent().parent().show(defaultSpeed);
	    $("#id_3-method_notes").parent().parent().show(defaultSpeed);
	} else {
	    $("#id_3-is_coordinate").prop("checked", false);
	    $("#id_3-peak_intensity").parent().parent().hide(defaultSpeed);
	    $("#id_3-peak_calling_method").parent().parent().hide(defaultSpeed);
	    $("#id_3-assay_conditions").parent().parent().hide(defaultSpeed);
	    $("#id_3-method_notes").parent().parent().hide(defaultSpeed);
	}
    });

    */
    
});


