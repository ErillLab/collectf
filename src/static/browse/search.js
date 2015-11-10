$(document).ready(function() {

    // experimental technique category change
    $('.category').on('change', function() {
        var id = $(this).prop('id');
        var selectValue = $(this).val();
        $('.'+id).css("display", "none");
        $('.'+id+'.'+selectValue).css("display", "block");
    });

    // unselect any selected techniques in not-selected category.
    $('input[type=submit]').on('click', function() {
        var selectValue1 = $('#cat1').val();
        var selectValue2 = $('#cat2').val();
        var selectValue3 = $('#cat3').val();
        $.each($('.cat1'), function() {
	        if (! $(this).hasClass(selectValue1) )
	            $(this).tree('unselectAll');
        });
        $.each($('.cat2'), function() {
	        if (! $(this).hasClass(selectValue2) )
	            $(this).tree('unselectAll');
        });
        $.each($('.cat3'), function() {
	        if (! $(this).hasClass(selectValue3) )
	            $(this).tree('unselectAll');
        });
    });

    // search tree
    $('.search_input').on('keyup', function() {
        var x = $(this).attr('id').replace('_search','');
        $('.'+x+':visible').tree('search', $(this).val());
        if($(this).val()=='')
	        $('.'+x+':visible').tree('collapseAll');
        return false;
    });

    // select all
    $('.select_box').change(function(ev) {
        var x = $(this).attr('id').replace('_select','');
        if ($(this).prop('checked') == true)
	        $('.'+x+':visible').tree('selectAll');
        else
	        $('.'+x+':visible').tree('unselectAll');

        return false;
    });


    // expand/collapse tree
    $('.toggle_link').on('click', function() {
        var x = $(this).attr('id').replace('_toggle','');
        $('.'+x+':visible').tree('toggleAll');
        return false;
    });


});
