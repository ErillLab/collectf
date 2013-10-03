// Sounds strange but this javascript is generated dynamically
//

{% load utils %}

var tf_data = {
    'child': [
	{% for TF_family in TF_families %}
	{
	    'text': '{{ TF_family.name }}',
	    'child': [
		{% for TF in TF_family.tf_set.all %}
		{
		    'text': '{{ TF.name }}',
                    'value': '{{ TF.TF_id }}'
		}{% if not forloop.end %},{% endif %}
		{% endfor %}
	    ]
	}{% if not forloop.last %},{% endif %}
	{% endfor %}
    ]
}

var species_data = {
    'child': [
	{% for phylum in phyla %}
	{
	    'text': '{{ phylum.name }}',
	    'child': [
		{% for class_ in phylum.taxonomy_set.all %}
		{
		    'text': '{{ class_.name }}',
		    'child': [
			{% for order in class_.taxonomy_set.all %}
			{
			    'text': '{{ order.name }}',
			    'child': [
				{% for family in order.taxonomy_set.all %}
				{
				    'text': '{{ family.name }}',
				    {% with genera=family.taxonomy_set.all %}
				    {% if genera %}
				    'child': [
					{% for genus in genera %}
					{
					    'text': '{{ genus.name }}',
					    {% with species_all=genus.taxonomy_set.all %}
					    {% if species_all %}
					    'child': [
						{% for species in species_all %}
						{
						    'text': '{{ species.name }}',
						    'value': '{{ species.pk }}'
						}{% if not forloop.last %},{% endif %}
						{% endfor %}
					    ]
					    {% else %}
					    'value': '{{ genus.pk }}'
					    {% endif %}
					    {% endwith %}
					}{% if not forloop.last %},{% endif %}
					{% endfor %}
				    ]
				    {% else %}
				    'value': '{{ family.pk }}'
				    {% endif %}
				    {% endwith %}
				}{% if not forloop.last %},{% endif %}
				{% endfor %}
			    ]
			}{% if not forloop.last %},{% endif %}
			{% endfor %}
		    ]
		}{% if not forloop.last %},{% endif %}
		{% endfor %}
	    ]
	}{% if not forloop.last %},{% endif %}
	{% endfor %}
    ]
}

var binding_techniques = {
    'child': [
	{% for subcategory in binding_techniques|get_keys %}
	{
	    'text': '{{ subcategory }}',
	    'child': [
		{% for technique in binding_techniques|get_item:subcategory %}
		{
		    'text': '{{ technique.name}}',
		    'value': '{{ technique.technique_id }}'
		}{% if not forloop.last %},{% endif %}
		{% endfor %}
	    ]
	}{% if not forloop.last %},{% endif %}
	{% endfor %}
    ]
}

var expression_techniques = {
    'child': [
	{% for subcategory in expression_techniques|get_keys %}
	{
	    'text': '{{ subcategory }}',
	    'child': [
		{% for technique in expression_techniques|get_item:subcategory %}
		{
		    'text': '{{ technique.name}}',
		    'value': '{{ technique.technique_id }}'
		}{% if not forloop.last %},{% endif %}
		{% endfor %}
	    ]
	}{% if not forloop.last %},{% endif %}
	{% endfor %}
    ]
}

var insilico_techniques = {
    'child': [
	{% for subcategory in insilico_techniques|get_keys %}
	{
	    'text': '{{ subcategory }}',
	    'child': [
		{% for technique in insilico_techniques|get_item:subcategory %}
		{
		    'text': '{{ technique.name}}',
		    'value': '{{ technique.technique_id }}'
		}{% if not forloop.last %},{% endif %}
		{% endfor %}
	    ]
	}{% if not forloop.last %},{% endif %}
	{% endfor %}
    ]
}

$('#tf_tree').tree({
    values: tf_data,
    type: 'checkbox',
    input: {name: 'tf_input'},
    toggleLI: true
});
$('#tf_tree').tree('collapseAll');
$('#species_tree').tree({
    values: species_data,
    type: 'checkbox',
    input: {name: 'species_input'},
    toggleLI: true
});
$('#species_tree').tree('collapseAll');

$('.binding').tree({
    values: binding_techniques,
    type: 'checkbox',
    toggleLI: true
});

$('.binding').tree('collapseAll');

$('.expression').tree({
    values: expression_techniques,
    type: 'checkbox',
    toggleLI: true
});

$('.expression').tree('collapseAll');

$('.insilico').tree({
    values: insilico_techniques,
    type: 'checkbox',
    toggleLI: true
});
$('.insilico').tree('collapseAll');

$('.category').on('change', function() {
    var id = $(this).prop('id');
    var selectValue = $(this).val();
    $('.'+id).css("display", "none");
    $('.'+id+'.'+selectValue).css("display", "block");
});

$('input[type=submit]').on('click', function() {
    var selectValue1 = $('#cat1').val();
    var selectValue2 = $('#cat2').val();
    var selectValue3 = $('#cat3').val();
    $.each($('.cat1'), function() {
	if (! $(this).hasClass(selectValue1) )
	    $(this).tree('unselectAll');
    });
    $.each($('.cat2'), function() {
	if (! $(this).hasClass(selectValue2) )
	    $(this).tree('unselectAll');
    });
    $.each($('.cat3'), function() {
	if (! $(this).hasClass(selectValue3) )
	    $(this).tree('unselectAll');
    });
})

$('.search_input').on('keyup', function() {
    //console.log($(this).val());
    var x = $(this).attr('id').replace('_search','');
    $('.'+x+':visible').tree('search', $(this).val());
    if($(this).val()=='')
	$('.'+x+':visible').tree('collapseAll');
    return false;
});

$('.select_box').change(function(ev) {
    var x = $(this).attr('id').replace('_select','');
    if ($(this).prop('checked') == true)
	$('.'+x+':visible').tree('selectAll');
    else
	$('.'+x+':visible').tree('unselectAll');

    return false;
});

$('.toggle_link').on('click', function() {
    var x = $(this).attr('id').replace('_toggle','');
    $('.'+x+':visible').tree('toggleAll');
    return false;
});
