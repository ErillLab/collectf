$(document).ready(function() {

    // onclick
    $("#id_1-TF_species_same").click(function(){
	$("#id_1-TF_species")[0].disabled = this.checked;
    });

    $("#id_1-site_species_same").click(function(){
	$("#id_1-site_species")[0].disabled = this.checked;
    });

    // when page is loaded, if checkboxes are checked
    $("#id_1-TF_species")[0].disabled = $("#id_1-TF_species_same")[0].checked;
    $("#id_1-site_species")[0].disabled = $("#id_1-site_species_same")[0].checked;
});