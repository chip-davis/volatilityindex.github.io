function youTubes_makeDynamic() {
    var $ytIframes = $('iframe[src*="youtube.com"]');
    $ytIframes.each(function (i,e) {
        var $ytFrame = $(e);
        var ytKey; var tmp = $ytFrame.attr('src').split(/\//); tmp = tmp[tmp.length - 1]; tmp = tmp.split('?'); ytKey = tmp[0];
        var $ytLoader = $('<div class="ytLoader">');
        $ytLoader.append($('<img class="cover" src="https://i.ytimg.com/vi/'+ytKey+'/hqdefault.jpg" width="350" height="222" style="vertical-align:top; position: relative; top: 0%;transform: translateY(-12%);" > '));
        $ytLoader.append($('<img class="playBtn" src="img/play_button.png">'));
        $ytLoader.data('$ytFrame',$ytFrame);
        $ytFrame.replaceWith($ytLoader);
        $ytLoader.click(function () {
            var $ytFrame = $ytLoader.data('$ytFrame');
            $ytFrame.attr('src',$ytFrame.attr('src')+'?autoplay=1');
            $ytLoader.replaceWith($ytFrame);
        });
    });
};
$(document).ready(function () {

    youTubes_makeDynamic()
    

});

$('#playButton').click(function () {
    $('#carousel-slider').carousel('cycle');
    $('#pauseButton').show();
    $('#playButton').hide();
});
$('#pauseButton').click(function () {
    $('#carousel-slider').carousel('pause');
    $('#pauseButton').hide();
    $('#playButton').show();
});
$('#carousel-slider').carousel({
    interval:12000,
    pause: "false"
});

function validateForm() {
    if (isEmpty(document.getElementById('data_2').value.trim())) {
        alert('NAME is required!');
        return false;
    }
    if (isEmpty(document.getElementById('data_4').value.trim())) {
        alert('EMAIL is required!');
        return false;
    }
    if (!validateEmail(document.getElementById('data_4').value.trim())) {
        alert('EMAIL must be a valid email address!');
        return false;
    }
    return true;
}
function isEmpty(str) { return (str.length === 0 || !str.trim()); }
function validateEmail(email) {
    var re = /^([\w-]+(?:\.[\w-]+)*)@((?:[\w-]+\.)*\w[\w-]{0,66})\.([a-z]{2,15}(?:\.[a-z]{2})?)$/i;
    return isEmpty(email) || re.test(email);
}

function changeTooltip(obj){
    event.preventDefault();
    $(obj).attr('data-original-title', 'Copied');
    $(obj).tooltip('show');
    $(obj).attr('data-original-title', 'Click to Copy');
}