$(document).ready(function() {
// onclick
    if($("#id_1-TF_species_same").length > 0) {
	$("#id_1-TF_species_same").click(function(){
	    $("#id_1-TF_species")[0].disabled = this.checked;
	});

	$("#id_1-site_species_same").click(function(){
	    $("#id_1-site_species")[0].disabled = this.checked;
	});

	// when page is loaded, if checkboxes are checked
	$("#id_1-TF_species")[0].disabled = $("#id_1-TF_species_same")[0].checked;
	$("#id_1-site_species")[0].disabled = $("#id_1-site_species_same")[0].checked;
    }

    // twitter bootstrap popover
    $('body').popover({
	selector: '[data-toggle="popover"]',
	trigger: 'hover',
	placement: 'right',
	html: 'true',
    });

    // scroll page down a little bit
    window.addEventListener("hashchange", function() { scrollBy(0, -50); });
});