function textAreaUpload() {
    $('#id_3-sites').each(function () {
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
}

function enable_qval_format() {
    $("#id_3-quantitative_data_format").parent().parent().toggle();
}

$(document).ready(function (){
    textAreaUpload();

    // Change the ul style (hide bullet points on the list of experimental techniques)
    $('ul').css('list-style-type', 'none');

    // Make a clear distinction between regular sites and high-throughput information
    $('.form-group:has(#id_3-peaks)').before('<div class="page-header"><h3>High-throughput data</h3></div>');
});
