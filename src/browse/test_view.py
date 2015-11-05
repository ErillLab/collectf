from django.shortcuts import render

def test_view(request):
    """Test view function."""
    return render(request, 'test_view.html')

