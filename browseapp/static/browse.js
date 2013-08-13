
var block_message = '<h1>Please wait.<br/><i class="icon-spinner icon-spin"></i></h1>';
var block_css 

function getReportLinksTF(type, id) {
    // Given string type (TF_family|TF) and id of the object, perform ajax to
    // retrieve list of reports
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI({message: block_message});
    $.getJSON('/browse_TF_all_reports_json/'+type+'/'+id, function(data) {
 	// fill description field
	$('#intro').remove();
	$('#TF_description').remove();
	$('#report_table').remove()
	
	var desc = '<h4>' + data['TF_name'] + '</h4>';
	desc += '<p>' + data['description'] + '</p>';

	var TF_desc = $('<div id="TF_description" class="boxed"></div>');
	TF_desc.append(desc);

	var report_table = $('<div id="report_table" class="boxed"></div>');
	if (data['list'].length > 0) {
	    var all_link = $('<p><a href="' + data['view_all_report_link'] + '">View all.</a></p>');
	    // get data and create a table of it
	    var table = $('<table></table>').addClass('table table-condensed table-striped');
            table.append('<thead><th>TF</th><th>Species</th><th>Report</th></thead>');
	    table = table.append('<tbody></tbody>');
	    // iterate over array elms
	    $.each(data['list'], function(key,value) {
		table.append('<tr>' +
			     '<td>' + value['TF_name'] + '</td>' +
			     '<td>' + value['species_name'] + '</td>' + 
			     '<td>' + '<a class="block_before_load" href="' + value['report_link'] + '">view</a></td>' +
			     '</tr>');
			    });
	    report_table.append(all_link).append(table);
	} else {
	    report_table.append('<p>No data</p>');
	}
	$('#browse').append(TF_desc).append(report_table);
    });
}

function getReportLinksSpecies(id, name) {
    // given species id, return all curations of speices under the taxonomy level
    // defined by id
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI({message: block_message});
    $.getJSON('/browse_tax_all_reports_json/'+id, function(data) {
 	// fill description field
	$('#intro').remove();
	$('#description').remove();
	$('#report_table').remove();

	var desc = $('<div id="description" class="boxed"></div>');
	get_wiki(name);


	var report_table = $('<div id="report_table" class="boxed"></div>');
	if (data['list'].length > 0) {
	    var all_link = '<a class="block_before_load" href="' + data['view_all_report_link'] + '"> View all reports.</a>';
	    // get data and create a table of it
	    var table = $('<table></table>').addClass('table table-condensed table-striped');
            table.append('<thead><th>TF</th><th>Species</th><th>Report</th></thead>');
	    table = table.append('<tbody></tbody>');
	    // iterate over array elms
	    $.each(data['list'], function(key,value) {
		table.append('<tr>' +
			     '<td>' + value['TF_name'] + '<a class="block_before_load" href=http://www.ncbi.nlm.nih.gov/protein/?term=' + value['TF_accession'] + '> [' + value['TF_accession'] + ']</a></td>' +
			     '<td>' + '<a class="block_before_load" href="http://www.ncbi.nlm.nih.gov/taxonomy/?term=' + value['species_name'] + '"> ' + value['species_name'] + '</a></td>' + 
			     '<td>' + '<a class="block_before_load" href="' + value['report_link'] + '"> [view]</a></td>' + 
			     '</tr>');
	    });
	    report_table.append(all_link).append(table);
	} else {
	    report_table.append('<p>No data.</p>');
	}
	$('#browse').append(desc).append(report_table);
    });
}

function get_wiki(name) {
    // given title, return first paragraph of wikipedia article
    console.log(name);
    $.ajax({
	url: 'http://en.wikipedia.org/w/api.php',
	data: {
	    action: 'parse',
	    //titles: name,
	    page: name,
	    format: 'json',
	    prop: 'text',
	},
	dataType:'jsonp',
	success: function(data) {
	    var wikiHTML = data.parse.text["*"];
	    var wikiDOM = $("<document>"+wikiHTML+"</document>");
	    var first_p = $(wikiDOM.children('p').get(0));

	    first_p.find('sup').remove();
	    first_p.find('a').each(function() {
		$(this)
		    .attr('href', 'http://en.wikipedia.org'+$(this).attr('href'))
		    .attr('target','wikipedia');
	    });
	    
	    var desc = $('#description');
	    desc.append(first_p)
		.append("<a href='http://en.wikipedia.org/wiki/"+name+"' target='wikipedia'>Read more on Wikipedia</a>");
	}
    });
}

function getReportLinksTechniques(type, id) {
    // given technique type and id, return all reports that have that technique in
    // it.
    $(document).ajaxStop($.unblockUI); // unblock when ajax activity stops
    $.blockUI({message: block_message});
    $.getJSON('/browse_techniques_all_reports_json/'+type+'/'+id, function(data) {
 	// fill description field
	$('#intro').remove();
	$('#description').remove();
	$('#report_table').remove();

	var desc = $('<div id="description" class="boxed"></div>');
	desc.append('<h4>' + data['title'] + '</h4>')
	    .append('<p>' + data['description'] + '</p>')
	    .append('<p><a class="block_before_load" href="' + data['view_all_report_link'] + '"> View all reports.</a></p>');
	
	var report_table = $('<div id="report_table" class="boxed"></div>');
	// get data and create a table of it
	var table = $('<table></table>').addClass('table table-condensed');
	table.append('<thead><th>TF</th><th>Species</th><th>Report</th></thead>');
	// iterate over array elms
	$.each(data['list'], function(key,value) {
	    table.append('<tr>' +
			 '<td>' + value['TF_name'] + ' <a href="http://www.ncbi.nlm.nih.gov/protein/?term=' + value['TF_accession'] + '">[' + value['TF_accession'] + ']</a></td>' +
			 '<td>' + '<a href="http://www.ncbi.nlm.nih.gov/taxonomy/?term=' + value['species_name'] + '"> ' + value['species_name'] + '</a></td>' + 
			 '<td>' + '<a class="block_before_load" href="' + value['report_link'] + '"> [view]</a></td>' + 
			 '</tr>');
	});
	report_table.append(table);
	$('#browse').append(desc).append(report_table);
    });
}

$(document).ready(function() {
    $(document).on('click', 'a.block_before_load', function(ev) {
	ev.preventDefault();
	console.log($(this).attr('href'));
	$.blockUI({message: block_message});
	window.location = $(this).attr('href');
    });

});