// Used in left frame of 'Browse by TF/taxonomy/techniques' pages.


function getWiki(name) {
    // Given title, returns first paragraph of Wikipedia article.
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
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
            console.log(data);
            var wikiDOM = $("<document>" + data.parse.text['*'] + "</document>");
            var first_p = $(wikiDOM.children('p').get(0));

            first_p.find('sup').remove();
            first_p.find('a').each(function () {
                $(this)
                    .attr('href', 'http://en.wikipedia.org' + $(this).attr('href'))
                    .attr('target', 'wikipedia');
            });
            var desc = $('#result_info_description');
            desc.text('');
            desc.append(first_p)
                .append("<a href='http://en.wikipedia.org/wiki/" +
                        name +
                        "' target='wikipedia'>Read more on Wikipedia</a>");
        }
    });
}

// Performs AJAX call to retrieve list of motif reports by given TF family.
function getMotifReportsByTFFamily(objectID) {
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: '/browse/get_results_by_TF_family/' + objectID,
        success: function (data) {
            $("#browse").html(data);
        }
    });
}

// Performs AJAX call to retrieve list of motif reports by given TF.
function getMotifReportsByTF(objectID) {
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: "/browse/get_results_by_TF/" + objectID,
        success: function (data) {
            $("#browse").html(data);
        }
    });
}

// Retrieves list of motif report links by the given taxonomy ID.
function getMotifReportsByTaxonomy(objectId) {
    // retrieve list of reports
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: '/browse/get_results_by_taxonomy/' + objectId,
        success: function (data) {
            $("#browse").html(data);
            getWiki($("#browse h1").text());
        }
    });
}

// Retrieves the list of motif report links by the given category.
function getMotifReportsByTechniqueCategory(technique_function, objectId) {
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: ('/browse/get_results_by_technique_category/' +
              technique_function + '/' + objectId),
        success: function (data) {
            $("#browse").html(data);
        }
    });
}

// Retrieves the list of motif report links by binding/expression.
function getMotifReportsByTechniqueFunction(technique_function) {
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: '/browse/get_results_by_technique_function/' + technique_function,
        success: function (data) {
            $("#browse").html(data);
        }
    });
}


// Retrieves the list of motif report links by the given technique.
function getMotifReportsByTechnique(objectId) {
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: '/browse/get_results_by_technique/' + objectId,
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
        // Switch toggle-icon, if present.
        $(this).find('i').toggleClass('fa-toggle-off fa-toggle-on');
    });
});
