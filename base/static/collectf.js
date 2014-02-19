function overrideBlockUIDefaults() {
    // http://malsup.com/jquery/block/#options
    $.blockUI.defaults.message = '<h2>Please wait.<br/><i class="fa fa-spinner fa-spin"></i></h2>';
    $.blockUI.defaults.css.border = "3px";
}

$(document).ready(function() {

    // Bootstrap tooltip settings
    $('body').tooltip({
        selector: '[data-toggle="tooltip"]',
        trigger: 'hover',
        placement: 'top',
        html: 'true'
    });

    // Bootstrap popover settings
    $('body').popover({
        selector: '[data-toggle="popover"]',
        trigger: 'hover',
        placement: 'right',
        html: 'true'
    });

    // bootstrap & django hack
    // alias for alert-error class -> alert-danger
    $('.alert-error').addClass('alert-danger');

    // override blockUI defaults
    overrideBlockUIDefaults();

})
