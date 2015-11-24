from django.shortcuts import render


def test_view(request):
    """Test view function."""
    return render(request, 'test_view.html')


def view_404(request):
    """Tests the 404 page."""
    return render(request, '404.html')


def view_500(request):
    """Tests the 500 page."""
    return render(request, '500.html')
