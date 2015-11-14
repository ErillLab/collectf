$(document).ready(function() {
    // Alias for alert-error class: alert-danger
    $('.alert-error').addClass('alert-danger');

    // Override blockUI defaults
    overrideBlockUIDefaults();

    // Activate dropdowns
    $('.dropdown-toggle').dropdown();
});

function overrideBlockUIDefaults() {
    // http://malsup.com/jquery/block/#options
    $.blockUI.defaults.message = '<h2>Please wait.<br/><i class="fa fa-spinner fa-spin"></i></h2>';
    $.blockUI.defaults.css.border = "3px";
}
