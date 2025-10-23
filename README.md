# アドバンストCG / Advanced CG
筑波大学情報学群情報メディア創成学類・情報科学類の授業「アドバンストCG」の第5回～第8回レポートのためのベースとなるコード公開用リポジトリです．
第1回～第4回については金森先生の[リポジトリ](https://github.com/yshhrknmr/AdvancedCG)を参照sてください．
This is a git repository for the lecture "Advanced Computer Graphics" at the University of Tsukuba, School of Informatics, College of Media Arts, Science and Technology. It includes the base code written in C++ for the reports.

## フォルダ構成など
- Advanced05～Advanced08: 第5回～第8回の授業課題のベースプログラムが置いてあります．中にはソースコードとVisual Studio用のプロジェクトファイルなどが入っています．Visual Studioから開く場合はこのフォルダ内のプロジェクトファイルではなく，リポジトリ直下に置いてあるAdvanced05Assignment_vs2022.sln (自分の環境などでVS2019を使っている場合は_vs2019.slnの方)を開いてください．
- Resources: 実行時フォルダで，実行時に参照するモデルファイルなどもここに置いてあります．
- include: 第5回～第8回の授業課題で必要となるヘッダファイルが置いてあります．Visual Studio 以外でコンパイルする際はこのディレクトリ (include) をインクルードフォルダとして追加してください．
- vc-x64-lib-common: Visual Studio 共通のライブラリファイル
- vc-x64-libs: Visual Studio のバージョン別のライブラリファイル

## トラブルシューティング
 - Visual Studioはgithub連携機能がありますが，このリポジトリには複数のVisual Studioソリューションとソースコード以外のファイルが多く含まれているのでうまくいかない可能性があります．別途リポジトリをクローンするかDownload ZIPで全体をダウンロード・解凍後，srcフォルダの各フォルダ内に含まれる*.slnファイルをVisual Studioで開くようにしてください．

 - Visual Studioで「error MSB8036: Windows SDK バージョン 10.0.xxxxxx.0 が見つかりませんでした。」というエラーが出てビルドできない．  
　⇒ 「プロジェクト」メニューから「プロジェクトの再ターゲット」で「Windows SDK バージョン:」 のところにその環境で対応するバージョンが出るので，問題なければそのまま「OK」をクリック
 - Visual Studioで実行しても特定のプロジェクトしか実行されない場合は，「プロジェクト」メニューから「スタートアップ プロジェクトの構成」で「現在の選択」にチェックしてください．

 - 自身のPC環境(Windows)で動かしたい人向けに講義資料フォルダ(SharePointの方)に「Visual Studio Community インストール方法」というPDFを用意してあるので参考にしてください．
