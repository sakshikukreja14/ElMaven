#include "numeric_treewidgetitem.h"
#include "PeakGroup.h"

static PeakGroup::ClassifiedLabel labelForString(const QString& labelString)
{
    if (labelString == "PeakGroup::ClassifiedLabel::Signal") {
        return PeakGroup::ClassifiedLabel::Signal;
    } else if (labelString == "PeakGroup::ClassifiedLabel::Noise") {
        return PeakGroup::ClassifiedLabel::Noise;
    } else if (labelString == "PeakGroup::ClassifiedLabel::Correlation") {
        return PeakGroup::ClassifiedLabel::Correlation;
    } else if (labelString == "PeakGroup::ClassifiedLabel::Pattern") {
        return PeakGroup::ClassifiedLabel::Pattern;
    } else if (labelString == "PeakGroup::ClassifiedLabel::CorrelationAndPattern") {
        return PeakGroup::ClassifiedLabel::CorrelationAndPattern;
    }
    return PeakGroup::ClassifiedLabel::None;
}

bool NumericTreeWidgetItem::operator<( const QTreeWidgetItem & other ) const{
    int sortCol = treeWidget()->sortColumn();

    // takes care of sorting based on PeakML class labels
    QVariant thisUserData = this->data(sortCol, Qt::UserRole);
    QVariant otherUserData = other.data(sortCol, Qt::UserRole);
    auto thisLabel = labelForString(thisUserData.value<QString>());
    auto otherLabel = labelForString(otherUserData.value<QString>());
    if (thisLabel != PeakGroup::ClassifiedLabel::None
        || otherLabel != PeakGroup::ClassifiedLabel::None) {
        return thisLabel < otherLabel;
    }

    QString thisText = text(sortCol);
    QString otherText = other.text(sortCol);

    QCollator collator;
    collator.setNumericMode(true);

    QRegExp exp("^[-+]?([0-9]*\\.?[0-9]+)[eE][-+]?([0-9]+)$");
    if (thisText.contains(exp)) {
        double mantissa = exp.cap(1).toDouble();
        double exponent = exp.cap(2).toDouble();
        double trueValue = mantissa * (pow(10, exponent));
        thisText = QString::fromStdString(to_string(trueValue).c_str());
    }
    if (otherText.contains(exp)) {
        double mantissa = exp.cap(1).toDouble();
        double exponent = exp.cap(2).toDouble();
        double trueValue = mantissa * (pow(10, exponent));
        otherText = QString::fromStdString(to_string(trueValue).c_str());
    }

    return collator.compare(thisText , otherText) < 0;
}
