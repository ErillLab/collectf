// Used in left frame of 'Browse by TF/taxonomy/techniques' pages.


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

function getMotifReportsByTFFamily(objectID) {
    $.ajax({
        type: "GET",
        url: '/browse/get_results_by_TF_family/' + objectID,
        success: function (data) {
            $("#browse").html(data);
        }
    });
}

function getMotifReportsByTF(objectID) {
    $.ajax({
        type: "GET",
        url: "/browse/get_results_by_TF/" + objectID,
        success: function (data) {
            $("#browse").html(data);
        }
    });
}


$(document).ready(function () {
    "use strict";
    var expandIcon = "fa fa-angle-right";
    var collapseIcon = "fa fa-angle-down";
    
    $('.node > a').on('click', function () {
        $(this).parent().next().toggle(200);
        if ($(this).parent().hasClass('expandable')) {
            $(this).parent().children('i').toggleClass(collapseIcon);
            $(this).parent().children('i').toggleClass(expandIcon);
        }
    });
    $('.node').each(function () {
        if ($(this).hasClass('expandable')) {
            $(this).next().hide();
            $(this).children('i').addClass(expandIcon);
            
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
                $(this).children('i').toggleClass(collapseIcon);
                $(this).children('i').toggleClass(expandIcon);
            }
        });
    });
});
