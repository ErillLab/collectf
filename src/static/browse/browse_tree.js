// Used in left frame of 'Browse by TF/taxonomy/techniques' pages.

$(function () {
    "use strict";
    var expandClass = "fa fa-angle-right";
    var collapseClass = "fa fa-angle-down";
    $('.node > a').on('click', function () {
        $(this).parent().next().toggle(200);
        if ($(this).parent().hasClass('expandable')) {
            $(this).parent().children('i').toggleClass(collapseClass);
            $(this).parent().children('i').toggleClass(expandClass);
        }
    });
    $('.node').each(function () {
        if ($(this).hasClass('expandable')) {
            $(this).next().hide();
            $(this).children('i').addClass(expandClass);
            
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
                $(this).children('i').toggleClass(collapseClass);
                $(this).children('i').toggleClass(expandClass);
            }
        });
    });


    function get_wiki(name) {
        // Given title, returns first paragraph of Wikipedia article.
        $.ajax({
            url: 'http://en.wikipedia.org/w/api.php',
            data: {
                action: 'parse',
                page: name,
                format: 'json',
                prop: 'text'
            },
            dataType: 'jsonp',
            success: function (data) {
                var wikiDOM = $("<document>" + data.parse.text['*'] + "</document>");
                var first_p = $(wikiDOM.children('p').get(0));

                first_p.find('sup').remove();
                first_p.find('a').each(function () {
                    $(this)
                        .attr('href', 'http://en.wikipedia.org' + $(this).attr('href'))
                        .attr('target', 'wikipedia');
                });
                var desc = $('#description');
                desc.text('');
                desc.append('<h4>' + name + '</h4>');
                desc.append(first_p)
                    .append("<a href='http://en.wikipedia.org/wiki/" +
                            name +
                            "' target='wikipedia'>Read more on Wikipedia</a>");
            }
        });
    }

    
});
