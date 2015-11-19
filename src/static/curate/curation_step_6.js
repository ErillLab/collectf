
function selectAllSites() {
    var checkboxes = $("[id$='_site']");
    checkboxes.prop('checked', !checkboxes.prop('checked'));
}

function applyToSelectedSites(techniqueId) {
    // Given a technique ID, check checkboxes of selected sites for that technique.
    var checkboxes = $("[id$='_site']");
    for(i=0; i<checkboxes.length; i++) {
        if (checkboxes[i].checked) {
            var t = i + '_technique_' + techniqueId;
            var technique = $('#id_6-' + t );
            technique.prop('checked', 'checked');
        }
    }
}

function applyTFFunctionToSelectedSites() {
    // Apply TF function to the selected sites.
    var checkboxes = $("[id$='_site']");
    var v = $("#TF_function_selected").val();
    for(i=0; i<checkboxes.length; i++) {
        if (checkboxes[i].checked) {
            $("#id_6-" + i + '_TF_function').val(v);
        }
    }
}

function applyTFTypeToSelectedSites() {
    // Apply TF function to the selected sites.
    var checkboxes = $("[id$='_site']");
    var v = $("#TF_type_selected").val();
    for(i=0; i<checkboxes.length; i++) {
        if (checkboxes[i].checked) {
            $("#id_6-" + i + '_TF_type').val(v);
        }
    }
}

function diagramsOnHover() {
    // make diagrams to be visible on hover
    var seqs = $('.sequence');

    for (i=0; i<seqs.length; i++) {
        var seq = $(seqs[i]);
        seq.popover({title:'',
                     content:seq.next().append(seq.next().next()),
                     html:true,
                     trigger:'hover',
                     placement:'right'});
        seq.popover('show');
        seq.popover('hide');

    }
}

function clearAll(techniqueId) {
    var t = '_technique_' + techniqueId;
    $("[id$=" + "'" + t + "']").removeAttr('checked');
}

var lastChecked = null;

$(document).ready(function() {
    // remove class from all elements
    $('div').removeClass('form-group');
    // Add it back to the button
    $('.btn').parent().parent().addClass('form-group');
    // remove class form-control
    $('.form-control').removeClass('form-control');
    // Hide gene diagrams by default
    diagramsOnHover();
    // shift select multiple checkboxes
    var checkboxes = $('[id^=id_6-][id$=_site]');
    console.log(checkboxes);
    $(checkboxes).click(function(e) {
        console.log(lastChecked);
        if (!lastChecked) {
            lastChecked = this;
            return;
        }
        if (e.shiftKey) {
            var start = $(checkboxes).index(this);
            var end = $(checkboxes).index(lastChecked);
            console.log(start + ' ' + end);
            $(checkboxes.slice(Math.min(start, end),
                               Math.max(start, end)+1))
                .prop('checked', $(lastChecked).is(':checked'));
        }
        lastChecked = this;
    });

});
