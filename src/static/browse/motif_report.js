// Motif report page utility functions.

function hideCol(colClass) {
    // Hides columns with the given class.
    $('table .'+colClass).each(function() {
        $(this).toggle();
    });
    $('ul#hiddenCols').append(
        '<li id="' + colClass +
		'"><a href="javascript:;" onclick="showCol(\'' +
	    colClass + '\');">Show ' +
        colClass.substr(0, colClass.length-4) +
        ' column</a></li>');
}

function showCol(columnClass){
    // Shows columns with the given class.
    $('table .'+columnClass).each(function(index) {
        $(this).show();
    });
    $('li#'+columnClass).remove();
}

$(document).ready(function() {
    // Regulation diagram hide/show
    $("#regulation_diagram_switch").on('click', function() {
        $(".regulation_diagram").toggle();
        return false;
    });

    // Hide 'genome' and 'TF' columns in the table
    hideCol('genome_col');
    hideCol('TF_col');
    
});
    

