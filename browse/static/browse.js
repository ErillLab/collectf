/*global $,document,window*/
"use strict";

if (!window.console) {
    window.console = {};
}
if (!window.console.log) {
    window.console.log = function () { };
}

function get_wiki(name) {
    // given title, return first paragraph of wikipedia article
    $.ajax({
        url: 'http://en.wikipedia.org/w/api.php',
        data: {
            action: 'parse',
            //titles: name,
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
            //console.log($(first_p));
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

// Retrieves list of motif report links by the given taxonomy ID.
function getReportLinksByTaxonomy(objectId) {
    // retrieve list of reports
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: '/browse/get_results_tax/' + objectId,
        success: function (data) {
            $("#browse").html(data);
            get_wiki($("#browse #description h4").text());
        }
    });
}

// Retrieves the list of motif report links by the given category.
function getReportLinksByTechniqueCategory(technique_function, objectId) {
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: ('/browse/get_results_technique_category/' +
              technique_function + '/' + objectId),
        success: function (data) {
            $("#browse").html(data);
        }
    });
}

// Retrieves the list of motif report links by binding/expression.
function getReportLinksByTechniqueFunction(technique_function) {
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: '/browse/get_results_technique_all/' + technique_function,
        success: function (data) {
            $("#browse").html(data);
        }
    });
}


// Retrieves the list of motif report links by the given technique.
function getReportLinksByTechnique(objectId) {
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: '/browse/get_results_technique/' + objectId,
        success: function (data) {
            $("#browse").html(data);
        }
    });
}

// Performs AJAX call to retrieve list of motif reports by given TF family.
function getReportLinksByTFFamily(objectID) {
    $(document).ajaxStop($.unblockUI); // unblock when AJAX activity stops.
    $.blockUI();
    $.ajax({
        type: "GET",
        url: "/browse/get_results_TF_family/" + objectID,
        success: function(data) {
            $("#browse").html(data);
        }
    });
}

function getReportLinksByTF(objectID) {
    // Given TF id, performs ajax toretrieve list of reports.
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI();
    $.ajax({
        type: "GET",
        url: '/browse/get_results_TF/' + objectID,
        success: function (data) {
            $("#browse").html(data);
        }
    });
}

$(document).ready(function () {

});
