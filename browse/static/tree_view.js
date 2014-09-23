// used in left frame of 'browse by TF/taxonomy/techniques' pages.

$(function () {
    "use strict";
    $('.node > a').on('click', function () {
        $(this).parent().next().toggle(200);
        if ($(this).parent().hasClass('expandable')) {
            $(this).parent().children('i').toggleClass('fa fa-angle-down');
            $(this).parent().children('i').toggleClass('fa fa-angle-right');
        }
    });
    $('.node').each(function () {
        if ($(this).hasClass('expandable')) {
            $(this).next().hide();
            $(this).children('i').addClass('fa fa-angle-right');
            
        }
            $(this).children('i').css('width', '1em');
            $(this).children('i').css('display', 'inline-block');
            $(this).css('padding-left', '20px');
            $(this).parent().children('ul').css('padding-left', '20px');
    });
    $('#toggle_all').on('click', function() {
        $('.node').each(function() {
            $(this).next().toggle(200);
            if ($(this).hasClass('expandable')) {
                $(this).children('i').toggleClass('fa fa-angle-down');
                $(this).children('i').toggleClass('fa fa-angle-right');
            }
        });
    });
});
