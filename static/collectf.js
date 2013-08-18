/*global $, jQuery, document, FileReader*/
"use strict";

$(document).ready(function () {
    // onclick

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

    // for all textarea fields, provide a file upload option
    $('#id_3-sites, #id_3-chip_data_extra_field').each(function () {
        $('<input type="file">').insertAfter($(this))
            .change(function () {
                var textarea = $(this).prev(),
                    file = $(this).get(0).files[0],
                    reader = new FileReader();
                reader.readAsText(file);
                reader.onload = function () {
                    textarea.val(reader.result);
                };
                reader.onerror = function (evt) {
                    textarea.val("error reading file.");
                };
            });
    });

    // FORM STEP 1
    if ($("#id_1-TF_species_same").length > 0) {
        $("#id_1-TF_species_same").click(function () {
            $("#id_1-TF_species")[0].disabled = this.checked;
        });
        $("#id_1-site_species_same").click(function () {
            $("#id_1-site_species")[0].disabled = this.checked;
        });
        // when page is loaded, if checkboxes are checked
        $("#id_1-TF_species")[0].disabled = $("#id_1-TF_species_same")[0].checked;
        $("#id_1-site_species")[0].disabled = $("#id_1-site_species_same")[0].checked;
    }
    $('#id_1-genome_accession').attr('autocomplete', 'off');
    $('#id_1-genome_accession').typeahead({
        'source': function (query, process) {
            $.getJSON('/get_genomes', function (data) {
                var arr = $.map(data, function (n, i) {
                    return n['genome_accession'] + ' - ' + n['organism'];
                });
                return process(arr);
            });
        },
        'updater': function (item) {
            return item.split(' ')[0];
        },
        'items': 20
    });

    $('#id_1-TF_accession').attr('autocomplete', 'off');
    $('#id_1-TF_accession').typeahead({
        'source': function (query, process) {
            $.getJSON('/get_TF_instances', function (data) {
                var arr = $.map(data, function (n, i) {
                    return n['protein_accession'] + ' - ' + n['description'];
                });
                return process(arr);
            });
        },
        'updater': function (item) {
            return item.split(' ')[0];
        },
        'items': 20,
    });


    // FORM STEP 3
    // scroll page down a little bit

    var defaultSpeed = 200;

    // JS for reported sites step of curation form.
    if ($('#id_3-is_chip_data').not(':checked')) {
        $('#id_3-assay_conditions').parent().parent().hide(defaultSpeed);
        $('#id_3-chip_method_notes').parent().parent().hide(defaultSpeed);
        $('#id_3-chip_data_extra_field').parent().parent().hide(defaultSpeed);
    }

    if ($('#id_3-is_chip_data').is(':checked')) {
        $('#id_3-assay_conditions').parent().parent().show(defaultSpeed);
        $('#id_3-chip_method_notes').parent().parent().show(defaultSpeed);
        $('#id_3-has_quantitative_data').attr('checked', true);
        if ($('#id_3-is_motif_associated').is(':checked')) {
            $('#id_3-chip_data_extra_field').parent().parent().show(defaultSpeed);
        }
    }

    $("#id_3-is_chip_data").on('change', function () {
        if ($(this).is(":checked")) {
            $("#id_3-assay_conditions").parent().parent().show(defaultSpeed);
            $("#id_3-chip_method_notes").parent().parent().show(defaultSpeed);
            $('#id_3-has_quantitative_data').attr('checked', true);
            if ($('#id_3-is_motif_associated').is(':checked')) {
                $('#id_3-chip_data_extra_field').parent().parent().show(defaultSpeed);
            }
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

    $('#id_3-has_quantitative_data').on('change', function () {
        if ($(this).is(':checked')) {
            $('#id_3-quantitative_data_format').parent().parent().show(defaultSpeed);
        } else {
            $('#id_3-quantitative_data_format').parent().parent().hide(defaultSpeed);
        }
    });

    if ($('#id_3-is_motif_associated').not(':checked')) {
        $('#id_3-chip_data_extra_field').parent().parent().hide(defaultSpeed);
    }

    if ($('#id_3-is_motif_associated').is(':checked')) {
        if ($('#id_3-is_chip_data').is(':checked')) {
            $('#id_3-chip_data_extra_field').parent().parent().show(defaultSpeed);
        }
    }

    $('#id_3-is_motif_associated').on('change', function () {
        if ($(this).is(':checked')) {
            if ($('#id_3-is_chip_data').is(':checked')) {
                $('#id_3-chip_data_extra_field').parent().parent().show(defaultSpeed);
            }
        } else {
            $('#id_3-chip_data_extra_field').parent().parent().hide(defaultSpeed);
        }
    });

    // FORM STEP 4-5-7

    $('div')
        .filter(function () {
            return this.id.match('div_id_4-');
        })
        .find('.controls').addClass('boxed');

    $('div')
        .filter(function () {
            return this.id.match('div_id_5-');
        })
        .find('.controls').addClass('boxed');

    $('div')
        .filter(function () {
            return this.id.match('div_id_7-');
        })
        .find('.controls').addClass('boxed');


});


